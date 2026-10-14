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


def surname(name):
    """Reduce any feed's author spelling to a bare surname: "Wang N", "Walters,
    L. C.", "Paul Bunk" and "Thomas PG" all become one word."""
    name = re.sub(r"\s*\(.*?\)", "", clean(name))
    if "," in name:                       # "Walters, L. C." -> "Walters"
        name = name.split(",")[0]
    else:
        parts = [p for p in name.split() if p]
        if not parts:
            return ""
        last = parts[-1]
        # Trailing initials ("Wang N", "Thomas PG") leave the surname before them.
        if len(parts) > 1 and (last.isupper() or re.fullmatch(r"(?:[A-Z]\.?)+", last)):
            parts = parts[:-1]
        # Keep particles: "van der Merwe PA" -> "van der Merwe".
        i = len(parts) - 1
        while i > 0 and parts[i - 1][:1].islower():
            i -= 1
        name = " ".join(parts[i:])
    return name.strip(".,").strip()


def short_authors(names):
    names = [surname(n) for n in names]
    names = [n for n in names if n]
    if not names:
        return ""
    return f"{names[0]} et al." if len(names) > 1 else names[0]


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
def _pubmed_fetch_ids(ids):
    items = []
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
        time.sleep(0.4)  # NCBI allows 3 requests/s without an API key
    return [it for it in items if it["title"]]


def _pubmed_search(queries, window, exclude, retmax):
    ids = []
    for query in queries:
        with Entrez.esearch(db="pubmed", term=f"({query}) AND {window} {exclude}",
                            retmax=retmax) as h:
            ids.extend(Entrez.read(h).get("IdList", []))
        time.sleep(0.4)
    return list(dict.fromkeys(ids))


def _pubmed_setup(cfg, start, end):
    pcfg = cfg["pubmed"]
    Entrez.email = os.getenv(pcfg.get("email_env", "ENTREZ_EMAIL"), "")
    window = f'("{start:%Y/%m/%d}"[EDAT] : "{end:%Y/%m/%d}"[EDAT])'
    return pcfg, window, pcfg.get("exclude", "")


def fetch_pubmed(cfg, start, end):
    pcfg, window, exclude = _pubmed_setup(cfg, start, end)
    try:
        ids = _pubmed_search(pcfg["queries"].values(), window, exclude,
                             pcfg.get("max_results", 600))
    except Exception as e:
        return [], f"PubMed search failed ({e})"
    try:
        return _pubmed_fetch_ids(ids), None
    except Exception as e:
        return [], f"PubMed fetch failed ({e})"


def fetch_pubmed_terms(cfg, terms, start, end):
    """A second PubMed pass for free-text terms taken from this week's posts."""
    pcfg, window, exclude = _pubmed_setup(cfg, start, end)
    if not terms:
        return [], None
    phrases = " OR ".join(f'"{t}"[tiab]' if " " in t else f"{t}[tiab]" for t in terms)
    try:
        ids = _pubmed_search([phrases], window, exclude, pcfg.get("term_max_results", 150))
        return _pubmed_fetch_ids(ids), None
    except Exception as e:
        return [], f"PubMed term search failed ({e})"


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


def _post_record(post, origin, since):
    """Normalise one post, or None when it is older than the window."""
    record = post.get("record", {})
    created = record.get("createdAt", "")
    try:
        if datetime.fromisoformat(created.replace("Z", "+00:00")) < since:
            return None
    except ValueError:
        pass
    external = (post.get("embed") or {}).get("external") or {}
    handle = post.get("author", {}).get("handle", "")
    rkey = post.get("uri", "").rsplit("/", 1)[-1]
    return {
        "feed": origin,
        "author": handle,
        "text": clean(record.get("text", ""))[:500],
        "link_title": clean(external.get("title", "")),
        "link": external.get("uri", ""),
        "likes": post.get("likeCount", 0),
        "reposts": post.get("repostCount", 0),
        "url": f"https://bsky.app/profile/{handle}/post/{rkey}",
    }


def _bsky_pages(base, headers, method, params, limit, key="feed"):
    """Yield items from a paginated app.bsky endpoint until `limit` is reached."""
    cursor, fetched = None, 0
    while fetched < limit:
        query = dict(params, limit=min(100, limit - fetched))
        if cursor:
            query["cursor"] = cursor
        r = requests.get(f"{base}/{method}", params=query, headers=headers, timeout=30)
        r.raise_for_status()
        data = r.json()
        batch = data.get(key, [])
        fetched += len(batch)
        yield from batch
        cursor = data.get("cursor")
        if not batch or not cursor:
            return


def bsky_session(_cache={}):
    """(base_url, headers) for Bluesky, signed in when .env has credentials."""
    if not _cache:
        try:
            token = _bsky_token()
        except Exception:
            token = None
        _cache["base"] = "https://bsky.social/xrpc" if token else "https://public.api.bsky.app/xrpc"
        _cache["headers"] = {"Authorization": f"Bearer {token}"} if token else {}
    return _cache["base"], _cache["headers"]


def fetch_follow_set(handles, cap=2000):
    """Everyone followed by the given accounts, as {did: handle}. Borrowing a
    colleague's follows is the quickest way to a trustworthy people list."""
    base, headers = bsky_session()
    people, problems = {}, []
    for handle in handles:
        try:
            for who in _bsky_pages(base, headers, "app.bsky.graph.getFollows",
                                   {"actor": handle}, cap, key="follows"):
                did = who.get("did")
                if did:
                    people[did] = who.get("handle", did)
        except Exception as e:
            problems.append(f"Bluesky follows of {handle} failed ({e})")
    return people, "; ".join(problems) or None


def fetch_people_posts(people, since, per_author=30):
    """Recent posts from each person, skipping replies."""
    base, headers = bsky_session()
    posts, failed = [], 0
    for did, handle in people.items():
        try:
            for entry in _bsky_pages(base, headers, "app.bsky.feed.getAuthorFeed",
                                     {"actor": did, "filter": "posts_no_replies"},
                                     per_author):
                rec = _post_record(entry.get("post", {}), "people", since)
                if rec:
                    rec["followed"] = True
                    posts.append(rec)
        except Exception:
            failed += 1
    problem = f"Bluesky could not read {failed} of {len(people)} people" if failed else None
    return posts, problem


def fetch_bluesky(cfg, start, people=None):
    """Posts from the people you follow, from custom feeds, lists and keyword
    searches. Feeds come and go, so the people and searches carry the section
    when a feed goes quiet."""
    bcfg = cfg.get("bluesky", {})
    since = datetime.combine(start, datetime.min.time(), tzinfo=timezone.utc)
    base, headers = bsky_session()
    cap = bcfg.get("max_posts_per_feed", 300)

    posts, problems, empty = [], [], []

    def collect(origin, method, params, key="feed", wrapped=True):
        got = 0
        try:
            for entry in _bsky_pages(base, headers, method, params, cap, key):
                post = entry.get("post", {}) if wrapped else entry
                rec = _post_record(post, origin, since)
                if rec:
                    posts.append(rec)
                    got += 1
        except Exception as e:
            problems.append(f"Bluesky {origin} failed ({e})")
            return
        if got == 0:
            empty.append(origin)

    if people:
        got, problem = fetch_people_posts(people, since,
                                          bcfg.get("posts_per_person", 30))
        posts.extend(got)
        if problem:
            problems.append(problem)
        if not got:
            empty.append(f"{len(people)} people followed")

    for name, uri in (bcfg.get("feeds") or {}).items():
        collect(f"feed '{name}'", "app.bsky.feed.getFeed", {"feed": uri})
    for name, uri in (bcfg.get("lists") or {}).items():
        collect(f"list '{name}'", "app.bsky.feed.getListFeed", {"list": uri})
    for query in bcfg.get("searches") or []:
        collect(f"search '{query}'", "app.bsky.feed.searchPosts",
                {"q": query, "sort": "latest", "since": since.strftime("%Y-%m-%dT%H:%M:%SZ")},
                key="posts", wrapped=False)

    if empty and not posts:
        problems.append("Bluesky returned no posts at all (" + ", ".join(empty) + ")")
    elif empty:
        problems.append("Bluesky sources with no posts this week: " + ", ".join(empty))

    handles = {h.lower() for h in (people or {}).values()}
    unique = {}
    for p in posts:                            # a post can appear in several sources
        p["followed"] = p.get("followed") or p["author"].lower() in handles
        unique.setdefault(p["url"], p)
    return list(unique.values()), "; ".join(problems) or None


def fetch_everything(cfg, start, end, people=None):
    """Returns (papers, bluesky_posts, problems)."""
    papers, problems = [], []
    for fetch in (lambda: fetch_pubmed(cfg, start, end),
                  lambda: fetch_biorxiv(cfg, start, end),
                  lambda: fetch_journals(cfg, start)):
        got, problem = fetch()
        papers.extend(got)
        if problem:
            problems.append(problem)
    posts, problem = fetch_bluesky(cfg, start, people)
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
