import os
from dotenv import load_dotenv

load_dotenv()

default_config = {
    "Entrez_email": os.getenv("ENTREZ_EMAIL", ""),
    "BIORXIV_RSS_URL": "https://connect.biorxiv.org/biorxiv_xml.php?subject=immunology",
    "handle_bsky": os.getenv("HANDLE_BSKY", ""),
    "bsky_app_password": os.getenv("BSKY_APP_PASSWORD", ""),
    "feed_uri": "at://did:plc:g3xk2d3jg6dma5glgwilyiyp/app.bsky.feed.generator/aaaknfoyzxnfa",
    "openai_api_key": os.getenv("OPENAI_API_KEY", ""),
    "use_pubmed": True,
    "use_biorxiv": True,
    "use_bluesky": True,
    "GPT_summary": False,
    "broad_terms": ["e.g.", "multiple sclerosis"
    ],
    "narrow_terms": ["e.g.", "T cells"
    ],
    "DAYS_BACK": 7,
    "FILTER_SCOPE": "title",
    "prompt1": (
        "You are an immunologist summarizing recent publications and posts. "
        "Focus on extracting key findings, trends, or notable publications relevant "
        "to T cells in multiple sclerosis and autoimmunity. Structure your summary "
        "the following way: write 2 sentences on highlights followed by a bullet point "
        "section for (if present) PubMed, biorxiv and Bluesky. If something is not relevant "
        "to the topic, ignore it. Be succinct and professional"
    )
}
