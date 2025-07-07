import re
import feedparser
import requests
import openai
import time
from config_file import default_config
from Bio import Entrez
from Bio import Medline
from datetime import datetime, timedelta

config = default_config.copy()

# -------------- REGEX ----------------
def build_patterns(config):
    broad_terms = config["broad_terms"]
    narrow_terms = config["narrow_terms"]
    broad_pattern = re.compile(r"\b(" + "|".join(broad_terms) + r")\b", re.IGNORECASE)
    narrow_pattern = re.compile(r"\b(" + "|".join(narrow_terms) + r")\b", re.IGNORECASE)
    return broad_pattern, narrow_pattern


# --------- FUNCTION DEFINITIONS ---------
def build_pubmed_query(terms, start_date, end_date, field="Title/Abstract"):
    formatted_terms = [
        f'"{term}"[{field}]' if " " in term else f'{term}[{field}]'
        for term in terms
    ]
    terms_query = "(" + " OR ".join(formatted_terms) + ")"
    date_filter = f'("{start_date}"[PDAT] : "{end_date}"[PDAT])'
    return f"{terms_query} AND {date_filter}"


def fetch_pubmed_results(config, max_retries=3, retry_delay=5):
    print("\n🔬 Fetching PubMed articles...")

    Entrez.email = config.get("Entrez_email", "")
    max_results = 50
    start_date = (datetime.now() - timedelta(days=config.get("DAYS_BACK", 30))).strftime("%Y/%m/%d")
    end_date = datetime.now().strftime("%Y/%m/%d")

    search_query = " OR ".join(config.get("broad_terms", []))
    if not search_query:
        print("⚠️ No broad terms defined. Skipping PubMed search.")
        return [], "⚠️ No broad terms defined.", 0

    print(f"🔎 PubMed query: '{search_query}', from {start_date} to {end_date}")

    attempt = 0
    while attempt < max_retries:
        try:
            handle = Entrez.esearch(
                db="pubmed",
                term=search_query,
                retmax=max_results,
                mindate=start_date,
                maxdate=end_date,
                datetype="pdat"
            )
            record = Entrez.read(handle)
            handle.close()
            break
        except Exception as e:
            attempt += 1
            print(f"⚠️ PubMed API error: {e} (Attempt {attempt}/{max_retries})")
            if attempt < max_retries:
                time.sleep(retry_delay)
            else:
                print("PubMed search failed after retries.")
                return [], "PubMed search failed after multiple attempts.", 0

    id_list = record.get("IdList", [])
    if not id_list:
        print("🔍 No PubMed results found.")
        return [], "No PubMed results found.", 0

    results = []
    try:
        fetch_handle = Entrez.efetch(db="pubmed", id=",".join(id_list), rettype="medline", retmode="text")
        records = Medline.parse(fetch_handle)

        for rec in records:
            title = rec.get("TI", "No Title")
            authors = ", ".join(rec.get("AU", []))
            published = rec.get("DP", "Unknown")
            doi = rec.get("AID", [""])[0].split(" ")[0] if "AID" in rec else ""
            abstract = rec.get("AB", "")

            if config.get("narrow_terms"):
                if not any(term.lower() in (title + abstract).lower() for term in config["narrow_terms"]):
                    continue

            results.append({
                "source": "PubMed",
                "title": title,
                "authors": authors,
                "published": published,
                "doi": doi,
                "abstract": abstract
            })

        fetch_handle.close()
        count = len(results)
        print(f"✅ Retrieved {count} PubMed results.")
        return results, "", count

    except Exception as e:
        print(f"Failed to fetch or parse PubMed details: {e}")
        return [], "Failed to fetch or parse PubMed details.", 0


# -------- BIORXIV PIPELINE --------
def parse_biorxiv_date(entry):
    date_str = None
    if hasattr(entry, 'dc_date'):
        date_str = entry.dc_date
    elif hasattr(entry, 'prism_publicationdate'):
        date_str = entry.prism_publicationdate

    if date_str:
        try:
            return datetime.strptime(date_str, "%Y-%m-%d")
        except Exception as e:
            print(f"Warning: could not parse date {date_str}: {e}")

    return None


def fetch_biorxiv_results(config):
    print("\n🧬 Fetching recent bioRxiv articles...")
    feed = feedparser.parse(config["BIORXIV_RSS_URL"])
    results = []

    if not feed.entries:
        print("No recent bioRxiv articles found.")
        return results, 0

    broad_pattern, narrow_pattern = build_patterns(config)

    for entry in feed.entries:
        published = parse_biorxiv_date(entry)
        if not published:
            print("⚠️ Skipping entry with no valid publication date.")
            continue

        title = entry.title
        abstract = entry.summary if hasattr(entry, "summary") else "No abstract available"

        authors = "Unknown"
        if hasattr(entry, "authors") and entry.authors:
            try:
                authors = ", ".join(a['name'] for a in entry.authors if 'name' in a)
            except Exception:
                authors = str(entry.authors)

        doi_str = "No DOI available"
        if hasattr(entry, "dc_identifier"):
            doi_str = entry.dc_identifier.replace("doi:", "").strip()

        text_to_search = title if config["FILTER_SCOPE"] == "title" else f"{title} {abstract}"

        if broad_pattern.search(text_to_search) and narrow_pattern.search(text_to_search):
            results.append({
                "source": "bioRxiv",
                "title": title,
                "authors": authors,
                "published": published.strftime("%Y-%m-%d"),
                "doi": doi_str,
                "abstract": abstract
            })

    count = len(results)
    if not results:
        print("No articles matched both broad and narrow filters.")

    return results, count


# -------- Bluesky PIPELINE --------
def get_session_token(config):
    url = "https://bsky.social/xrpc/com.atproto.server.createSession"
    payload = {"identifier": config["handle_bsky"], "password": config["bsky_app_password"]}
    resp = requests.post(url, json=payload)
    resp.raise_for_status()
    return resp.json()["accessJwt"]


def fetch_custom_bluesky_feed(config, limit=20):
    print("\n💬 Fetching Bluesky posts...")
    token = get_session_token(config)
    url = "https://bsky.social/xrpc/app.bsky.feed.getFeed"
    headers = {"Authorization": f"Bearer {token}"}
    params = {"feed": config["feed_uri"], "limit": limit}

    resp = requests.get(url, headers=headers, params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    posts = data.get("feed", [])
    if not posts:
        print("No posts found in the custom feed.")
        return results, 0

    broad_pattern, narrow_pattern = build_patterns(config)

    for post in posts:
        record = post.get("post", {}).get("record", {})
        text = record.get("text", "")

        if broad_pattern.search(text) and narrow_pattern.search(text):
            uri = post.get("post", {}).get("uri", "")
            created_at = record.get("createdAt")
            author_handle = post.get("post", {}).get("author", {}).get("handle", "unknown")

            published = "unknown"
            if created_at:
                try:
                    dt = datetime.fromisoformat(created_at.rstrip("Z"))
                    published = dt.strftime("%Y-%m-%d %H:%M:%S")
                except Exception:
                    published = created_at

            results.append({
                "source": "Bluesky",
                "title": f"Bluesky Post by @{author_handle}",
                "authors": author_handle,
                "published": published,
                "doi": uri,
                "abstract": text
            })

    count = len(results)
    if not results:
        print("No posts matched both broad and narrow filters.")

    return results, count


def summarize_results(all_results, config):
    summary_input = ""
    for entry in all_results:
        summary_input += (
            f"Source: {entry.get('source', 'Unknown')}\n"
            f"Title: {entry.get('title', 'N/A')}\n"
            f"Authors: {entry.get('authors', 'N/A')}\n"
            f"Published: {entry.get('published', 'N/A')}\n"
            f"DOI/Link: {entry.get('doi', 'N/A')}\n"
            f"Content: {entry.get('abstract', 'N/A')}\n"
            f"{'-'*50}\n"
        )

    openai.api_key = config["openai_api_key"]

    response = openai.chat.completions.create(
        model="gpt-4.1-mini",
        messages=[
            {"role": "system", "content": config["prompt1"]},
            {"role": "user", "content": f"Summarize the following literature and posts:\n\n{summary_input}"}
        ],
        max_tokens=800,
        temperature=0.3
    )

    return response.choices[0].message.content


# -------- New unified fetch function --------
def fetch_all_results(config):
    all_results = []
    grouped_results = {"PubMed": [], "bioRxiv": [], "Bluesky": []}
    searched_sources = []

    # PubMed
    if config.get("use_pubmed", False):
        searched_sources.append("PubMed")
        pubmed_results, pubmed_msg, pubmed_count = fetch_pubmed_results(config)
        grouped_results["PubMed"] = pubmed_results
    else:
        pubmed_count = 0

    # bioRxiv
    if config.get("use_biorxiv", False):
        searched_sources.append("bioRxiv")
        biorxiv_results, biorxiv_count = fetch_biorxiv_results(config)
        grouped_results["bioRxiv"] = biorxiv_results
    else:
        biorxiv_count = 0

    # Bluesky
    if config.get("use_bluesky", False):
        searched_sources.append("Bluesky")
        bluesky_results, bluesky_count = fetch_custom_bluesky_feed(config)
        grouped_results["Bluesky"] = bluesky_results
    else:
        bluesky_count = 0

    for src in searched_sources:
        all_results.extend(grouped_results[src])

    counts = {
        "PubMed": pubmed_count,
        "bioRxiv": biorxiv_count,
        "Bluesky": bluesky_count,
    }

    return all_results, grouped_results, counts, searched_sources
