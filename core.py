import re
import feedparser
import requests
import time
import gc
import os
from config_file import default_config
from Bio import Entrez
from Bio import Medline
from datetime import datetime, timedelta
import openai
from openai import OpenAIError, OpenAI
from math import ceil
from typing import List, Dict, Generator, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

config = default_config.copy()
client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))


# -------------- REGEX ----------------
def build_patterns(config):
    broad_terms = config["broad_terms"]
    narrow_terms = config["narrow_terms"]
    broad_pattern = re.compile(r"\b(" + "|".join(broad_terms) + r")\b", re.IGNORECASE)
    narrow_pattern = re.compile(r"\b(" + "|".join(narrow_terms) + r")\b", re.IGNORECASE)
    return broad_pattern, narrow_pattern


# --------- FUNCTION DEFINITIONS ---------
def build_pubmed_boolean_query(broad_terms, narrow_terms, start_date, end_date, field="Title/Abstract"):
    # Map config values to proper PubMed field tags
    field_map = {
        "title": "Title",
        "title/abstract": "Title/Abstract"
    }
    pubmed_field = field_map.get(field.lower(), "Title/Abstract")

    def format_terms(terms):
        return " OR ".join(
            f'"{term}"[{pubmed_field}]' if " " in term else f'{term}[{pubmed_field}]'
            for term in terms
        )

    broad_query = f"({format_terms(broad_terms)})" if broad_terms else ""
    narrow_query = f"({format_terms(narrow_terms)})" if narrow_terms else ""

    if broad_query and narrow_query:
        terms_query = f"{broad_query} AND {narrow_query}"
    elif broad_query:
        terms_query = broad_query
    elif narrow_query:
        terms_query = narrow_query
    else:
        terms_query = ""

    date_filter = f'("{start_date}"[PDAT] : "{end_date}"[PDAT])'
    if terms_query:
        return f"{terms_query} AND {date_filter}"
    else:
        return date_filter


def fetch_pubmed_results(config, max_retries=3, retry_delay=5):
    print("\n🔬 Fetching PubMed articles...")

    Entrez.email = config.get("Entrez_email", "")
    max_results = 500
    source_name = "PubMed"
    fetch_type = "medline"
    fetch_mode = "text"

    # Dates
    days_back = config.get("DAYS_BACK", 30)
    start_date = (datetime.now() - timedelta(days=days_back)).strftime("%Y/%m/%d")
    end_date = datetime.now().strftime("%Y/%m/%d")

    # Terms and field
    broad_terms = config.get("broad_terms", [])
    narrow_terms = config.get("narrow_terms", [])
    field = config.get("FILTER_SCOPE", "title")  # default is "title"

    # Build the query
    search_query = build_pubmed_boolean_query(broad_terms, narrow_terms, start_date, end_date, field=field)
    if not search_query:
        warning = "⚠️ No search terms defined. Skipping PubMed search."
        print(warning)
        return [], warning, 0

    print(f"🔎 PubMed query: '{search_query}'")

    # Search
    for attempt in range(1, max_retries + 1):
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
            print(f"⚠️ PubMed API error: {e} (Attempt {attempt}/{max_retries})")
            if attempt < max_retries:
                time.sleep(retry_delay)
            else:
                return [], "❌ PubMed search failed after multiple attempts.", 0

    id_list = record.get("IdList", [])
    if not id_list:
        print("🔍 No PubMed results found.")
        return [], "No PubMed results found.", 0

    # Fetch and parse details
    results = []
    try:
        with Entrez.efetch(db="pubmed", id=",".join(id_list), rettype=fetch_type, retmode=fetch_mode) as fetch_handle:
            for rec in Medline.parse(fetch_handle):
                title = rec.get("TI", "No Title")
                authors = ", ".join(rec.get("AU", []))
                published = rec.get("DP", "Unknown")
                abstract = rec.get("AB", "")
                doi = ""

                for aid in rec.get("AID", []):
                    if "doi" in aid.lower():
                        doi = "https://doi.org/" + aid.split(" ")[0]
                        break

                results.append({
                    "source": source_name,
                    "title": title,
                    "authors": authors,
                    "published": published,
                    "doi": doi,
                    "abstract": abstract
                })

        count = len(results)
        print(f"✅ Retrieved {count} PubMed results.")
        return results, "", count

    except Exception as e:
        print(f"❌ Failed to fetch or parse PubMed details: {e}")
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
    subject = config.get("biorxiv_subject", "Immunology")
    rss_url = f"https://connect.biorxiv.org/biorxiv_xml.php?subject={subject}"
    feed = feedparser.parse(rss_url)

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


def fetch_custom_bluesky_feed(config, limit=100):
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


# -------- GPT PIPELINE --------
def batch_generator(items: List[dict], batch_size: int) -> Generator[List[dict], None, None]:
    """
    Yield successive batches of entries from the list.
    """
    for i in range(0, len(items), batch_size):
        yield items[i: i + batch_size]


def format_entries(entries: List[dict]) -> str:
    """
    Format a list of entries into a plain text block for OpenAI prompt.
    """
    formatted_chunks = []
    for entry in entries:
        chunk = "\n".join([
            f"Title: {entry.get('title', 'N/A')}",
            f"Authors: {entry.get('authors', 'N/A')}",
            f"Published: {entry.get('published', 'N/A')}",
            f"Content: {entry.get('abstract', 'N/A')}",
            "-" * 10,
        ])
        formatted_chunks.append(chunk)
    return "\n".join(formatted_chunks)


def call_openai(prompt_text: str, api_key: str, system_prompt: str = "You are a helpful assistant.",
                timeout: int = 120, max_retries: int = 3, retry_delay: float = 2.0) -> str:
    openai.api_key = api_key
    for attempt in range(1, max_retries + 1):
        try:
            response = openai.chat.completions.create(
                model="gpt-4.1-mini",
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": prompt_text},
                ],
                max_tokens=800,
                temperature=0.3,
                timeout=timeout,
                stop=["[END SUMMARY]"],
            )
            return response.choices[0].message.content.strip()

        except OpenAIError as e:
            logger.warning(f"OpenAI API call failed on attempt {attempt}/{max_retries}: {e}")
            if attempt == max_retries:
                logger.error("Max retries reached. Returning empty string.")
                return ""
            time.sleep(retry_delay * attempt)

        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            return ""


def generate_terms_from_interest(interest_sentence, config):
    prompt = (
        "This is a sentence that describes the scientist's interests:\n"
        f"\"{interest_sentence}\"\n\n"
        "Generate up to 10 highly specific search terms based on the following research interest. The terms will be "
        "used for a PubMed search. Strongly prefer one-word terms but do not fuse words. Only include terms that are precise and "
        "relevant while skipping vague or overly broad ones. Return a single comma-separated list with no numbering, "
        "explanations, or additional text. "
    )

    response = client.chat.completions.create(
        model="gpt-4.1-mini",
        messages=[
            {"role": "system", "content": "You are a biomedical literature expert."},
            {"role": "user", "content": prompt},
        ],
        temperature=0.7,
        max_tokens=300,
    )

    return response.choices[0].message.content


def summarize_batch(entries, source_name, prompt_intro, api_key, system_prompt):
    source_header = f"--- Source: {source_name} ---\n\n"
    prompt = prompt_intro + source_header + format_entries(entries)
    summary = call_openai(prompt, api_key, system_prompt)
    gc.collect()
    return summary


def summarize_results(
        all_results: List[dict],
        grouped_results: Dict[str, List[dict]],
        config: dict,
        batch_size: int = 5,
        max_workers: int = 3,
) -> str:
    """
    Summarize all results by batching large sources and combining summaries.
    Uses concurrent threads to speed up PubMed batch summarizations.
    Frees memory aggressively after summarization.
    """

    api_key = config["openai_api_key"]
    total_results = len(all_results)
    final_summaries = []
    generic_system_prompt = "You are a scientist."

    batch_prompt_template = (
        "Summarize the following literature in bullet points, max 200 words. "
        "Include key findings, consensus, and cite sources explicitly as shown. "
        "Be concise and omit irrelevant information. "
        "End the summary with the phrase: [END SUMMARY]\n\n"
    )

    one_shot_prompt_template = batch_prompt_template

    combined_prompt_template = (
        "You are given summaries from different sources. "
        "Focus on overall trends and consensus, avoid repetition, and keep citations. "
        "Write max 200 words. "
        "End the summary with the phrase: [END SUMMARY]\n\n"
    )

    if total_results == 0:
        return "No literature results found to summarize."

    if total_results < 12:
        prompt = one_shot_prompt_template + format_entries(all_results)
        final_summaries.append(call_openai(prompt, api_key, generic_system_prompt))
        gc.collect()
    else:
        # PubMed batches (parallel)
        pubmed_entries = grouped_results.get("PubMed", [])
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for i, batch in enumerate(batch_generator(pubmed_entries, batch_size), start=1):
                futures.append(
                    executor.submit(
                        summarize_batch, batch, f"PubMed batch {i}", batch_prompt_template, api_key,
                        generic_system_prompt
                    )
                )
            for future in as_completed(futures):
                summary = future.result()
                final_summaries.append(summary)
        gc.collect()

        # bioRxiv (sequential)
        biorxiv_entries = grouped_results.get("bioRxiv", [])
        if biorxiv_entries:
            summary = summarize_batch(biorxiv_entries, "bioRxiv", one_shot_prompt_template, api_key,
                                      generic_system_prompt)
            final_summaries.append(f"--- bioRxiv summary ---\n{summary}")
            gc.collect()

        # Bluesky (sequential)
        bluesky_entries = grouped_results.get("Bluesky", [])
        if bluesky_entries:
            summary = summarize_batch(bluesky_entries, "Bluesky", one_shot_prompt_template, api_key,
                                      generic_system_prompt)
            final_summaries.append(f"--- Bluesky summary ---\n{summary}")
            gc.collect()

    # Combine summaries using user-defined prompt
    combined_text = "\n\n".join(final_summaries)

    # 🔧 Combine interest_sentence with prompt1, if interest sentence exists
    user_interest = config.get("interest_sentence", "").strip()
    user_prompt = config.get("prompt1", "").strip()
    if user_interest:
        combined_user_prompt = f"'{user_interest}'.\n\n{user_prompt}"
    else:
        combined_user_prompt = user_prompt

    final_prompt = combined_user_prompt + combined_text
    overall_summary = call_openai(final_prompt, api_key, combined_user_prompt)

    return overall_summary


# -------- Collection PIPELINE --------
def fetch_all_results(config: dict) -> Tuple[List[dict], Dict[str, List[dict]], Dict[str, int], List[str]]:
    """
    Fetch results from configured sources and group them.
    Returns:
        all_results: flat list of all entries
        grouped_results: dict with source keys and list of entries
        counts: dict of counts per source
        searched_sources: list of sources actually searched
    """
    all_results = []
    grouped_results = {"PubMed": [], "bioRxiv": [], "Bluesky": []}
    searched_sources = []

    if config.get("use_pubmed", False):
        searched_sources.append("PubMed")
        pubmed_results, _, pubmed_count = fetch_pubmed_results(config)
        grouped_results["PubMed"] = pubmed_results
    else:
        pubmed_count = 0

    if config.get("use_biorxiv", False):
        searched_sources.append("bioRxiv")
        biorxiv_results, biorxiv_count = fetch_biorxiv_results(config)
        grouped_results["bioRxiv"] = biorxiv_results
    else:
        biorxiv_count = 0

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
