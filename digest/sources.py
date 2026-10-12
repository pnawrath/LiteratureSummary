"""Fetchers for the weekly digest. Each returns (items, problem) where problem is
None on success or a short message that the issue prints instead of silently
dropping a source."""
import html
import os
import re
import time
from datetime import date, datetime, timedelta, timezone

import feedparser
import requests
from Bio import Entrez, Medline

USER_AGENT = "LiteratureDigest/1.0 (personal research digest)"
TAG_RE = re.compile(r"<[^>]+>")


def clean(text):
    return re.sub(r"\s+", " ", html.unescape(TAG_RE.sub(" ", text or ""))).strip()


def make_item(source, title, abstract="", authors="", published="", doi="", url="",
              journal="", kind="paper", general=False, pmid=""):
    doi = normalize_doi(doi)
    if not url and doi:
        url = f"https://doi.org/{doi}"
    return {
        "source": source,
        "journal": journal or source,
        "kind": kind,
        "general": general,
        "title": clean(title),
        "abstract": clean(abstract),
        "authors": authors,
        "published": published,
        "doi": doi,
        "pmid": pmid,
        "url": url,
    }


def normalize_doi(doi):
    doi = (doi or "").strip().lower()
    doi = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:)", "", doi)
    return doi if doi.startswith("10.") else ""


def short_authors(names):
    names = [n for n in names if n]
    if not names:
        return ""
    return names[0] if len(names) == 1 else f"{names[0]} et al."


def _get(url, **kwargs):
    last = None
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=60, headers={"User-Agent": USER_AGENT}, **kwargs)
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            last = e
            time.sleep(3 * (attempt + 1))
    raise last


# -------- PubMed --------
def fetch_pubmed(cfg, start, end):
    pcfg = cfg["pubmed"]
    Entrez.email = os.getenv(pcfg.get("email_env", "ENTREZ_EMAIL"), "")
    window = f'("{start:%Y/%m/%d}"[EDAT] : "{end:%Y/%m/%d}"[EDAT])'
    exclude = pcfg.get("exclude", "")
    ids = []
    try:
        for query in pcfg["queries"].values():
            term = f"({query}) AND {window} {exclude}"
            with Entrez.esearch(db="pubmed", term=term, retmax=pcfg.get("max_results", 600)) as h:
                ids.extend(Entrez.read(h).get("IdList", []))
            time.sleep(0.4)  # NCBI allows 3 requests/s without an API key
    except Exception as e:
        return [], f"PubMed search failed ({e})"

    ids = list(dict.fromkeys(ids))
    items = []
    try:
        for i in range(0, len(ids), 200):
            with Entrez.efetch(db="pubmed", id=",".join(ids[i:i + 200]),
                               rettype="medline", retmode="text") as h:
                for rec in Medline.parse(h):
                    doi = next((a.split(" ")[0] for a in rec.get("AID", []) if "[doi]" in a), "")
                    pmid = rec.get("PMID", "")
                    items.append(make_item(
                        "PubMed", rec.get("TI", ""), rec.get("AB", ""),
                        short_authors(rec.get("AU", [])), rec.get("DP", ""), doi,
                        url="" if doi else f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        journal=rec.get("TA", "PubMed"), pmid=pmid,
                    ))
            time.sleep(0.4)
    except Exception as e:
        return items, f"PubMed fetch stopped early ({e})"
    return [it for it in items if it["title"]], None


# -------- bioRxiv --------
def fetch_biorxiv(cfg, start, end):
    bcfg = cfg["biorxiv"]
    items, cursor = [], 0
    try:
        while True:
            url = bcfg["api_url"].format(start=f"{start:%Y-%m-%d}", end=f"{end:%Y-%m-%d}", cursor=cursor)
            data = _get(url).json()
            batch = data.get("collection", [])
            for rec in batch:
                # The API lists every version; keep the first version only.
                if str(rec.get("version", "1")) != "1":
                    continue
                authors = [a.strip() for a in rec.get("authors", "").split(";")]
                items.append(make_item(
                    "bioRxiv", rec.get("title", ""), rec.get("abstract", ""),
                    short_authors(authors), rec.get("date", ""), rec.get("doi", ""),
                    journal="bioRxiv", kind="preprint",
                ))
            total = int((data.get("messages") or [{}])[0].get("total", 0) or 0)
            cursor += len(batch)
            if not batch or cursor >= total:
                break
        if items:
            return items, None
        problem = "bioRxiv API returned nothing"
    except Exception as e:
        problem = f"bioRxiv API failed ({e})"

    # Fallback: RSS (only the latest ~30 preprints).
    feed_items, rss_problem = fetch_rss("bioRxiv", bcfg["rss_fallback"], start, kind="preprint")
    if feed_items:
        return feed_items, f"{problem}; used the RSS fallback, which only lists the latest preprints"
    return [], f"{problem}; RSS fallback: {rss_problem or 'no entries'}"


# -------- Journal tables of contents --------
def _entry_date(entry):
    for key in ("published_parsed", "updated_parsed"):
        if entry.get(key):
            return date(*entry[key][:3])
    for key in ("dc_date", "prism_publicationdate", "prism_coverdate"):
        raw = entry.get(key)
        if raw:
            try:
                return datetime.fromisoformat(raw[:10]).date()
            except ValueError:
                pass
    return None


def _entry_doi(entry):
    for key in ("prism_doi", "dc_identifier", "id", "link"):
        doi = normalize_doi(entry.get(key, ""))
        if doi:
            return doi
        m = re.search(r"(10\.\d{4,9}/[^\s?#]+)", entry.get(key, "") or "")
        if m:
            return normalize_doi(m.group(1))
    return ""


def fetch_rss(name, url, start, kind="paper", general=False):
    try:
        feed = feedparser.parse(_get(url).content)
    except Exception as e:
        return [], f"{name} feed failed ({e})"
    if not feed.entries:
        return [], f"{name} feed returned no entries (URL may have changed)"
    items = []
    for e in feed.entries:
        d = _entry_date(e)
        if d and d < start:
            continue
        authors = [a.get("name", "") for a in e.get("authors", [])] or [e.get("author", "")]
        items.append(make_item(
            name, e.get("title", ""), e.get("summary", "") or e.get("description", ""),
            short_authors(authors), d.isoformat() if d else "", _entry_doi(e),
            url=e.get("link", ""), journal=name, kind=kind, general=general,
        ))
    return [it for it in items if it["title"]], None


def fetch_journals(cfg, start):
    items, problems = [], []
    for j in cfg.get("journals", []):
        got, problem = fetch_rss(j["name"], j["url"], start, general=j.get("general", False))
        items.extend(got)
        if problem:
            problems.append(problem)
    return items, "; ".join(problems) or None


# -------- Bluesky --------
def _bsky_token():
    handle, password = os.getenv("BSKY_HANDLE"), os.getenv("BSKY_APP_PASSWORD")
    if not (handle and password):
        return None
    r = requests.post("https://bsky.social/xrpc/com.atproto.server.createSession",
                      json={"identifier": handle, "password": password}, timeout=30)
    r.raise_for_status()
    return r.json()["accessJwt"]


def fetch_bluesky(cfg, start):
    bcfg = cfg.get("bluesky", {})
    since = datetime.combine(start, datetime.min.time(), tzinfo=timezone.utc)
    try:
        token = _bsky_token()
    except Exception as e:
        return [], f"Bluesky login failed ({e})"
    base = "https://bsky.social/xrpc" if token else "https://public.api.bsky.app/xrpc"
    headers = {"Authorization": f"Bearer {token}"} if token else {}

    posts, problems = [], []
    for feed_name, uri in bcfg.get("feeds", {}).items():
        cursor, fetched = None, 0
        try:
            while fetched < bcfg.get("max_posts_per_feed", 300):
                params = {"feed": uri, "limit": 100}
                if cursor:
                    params["cursor"] = cursor
                r = requests.get(f"{base}/app.bsky.feed.getFeed", params=params,
                                 headers=headers, timeout=30)
                r.raise_for_status()
                data = r.json()
                batch = data.get("feed", [])
                fetched += len(batch)
                too_old = False
                for entry in batch:
                    post = entry.get("post", {})
                    record = post.get("record", {})
                    created = record.get("createdAt", "")
                    try:
                        if datetime.fromisoformat(created.replace("Z", "+00:00")) < since:
                            too_old = True
                            continue
                    except ValueError:
                        pass
                    external = (post.get("embed") or {}).get("external") or {}
                    handle = post.get("author", {}).get("handle", "")
                    rkey = post.get("uri", "").rsplit("/", 1)[-1]
                    posts.append({
                        "feed": feed_name,
                        "author": handle,
                        "text": clean(record.get("text", ""))[:500],
                        "link_title": clean(external.get("title", "")),
                        "link": external.get("uri", ""),
                        "likes": post.get("likeCount", 0),
                        "reposts": post.get("repostCount", 0),
                        "url": f"https://bsky.app/profile/{handle}/post/{rkey}",
                    })
                cursor = data.get("cursor")
                if not batch or not cursor or too_old:
                    break
        except Exception as e:
            problems.append(f"Bluesky feed '{feed_name}' failed ({e})")

    # Same post can appear in several feeds.
    unique = {p["url"]: p for p in posts}
    return list(unique.values()), "; ".join(problems) or None


def fetch_everything(cfg, start, end):
    """Returns (papers, bluesky_posts, problems)."""
    papers, problems = [], []
    for fetch in (lambda: fetch_pubmed(cfg, start, end),
                  lambda: fetch_biorxiv(cfg, start, end),
                  lambda: fetch_journals(cfg, start)):
        got, problem = fetch()
        papers.extend(got)
        if problem:
            problems.append(problem)
    posts, problem = fetch_bluesky(cfg, start)
    if problem:
        problems.append(problem)
    return papers, posts, problems


def window(cfg, last_run):
    end = date.today()
    if last_run:
        start = max(date.fromisoformat(last_run), end - timedelta(days=cfg["max_days_back"]))
    else:
        start = end - timedelta(days=cfg["days_back"])
    return start, end
