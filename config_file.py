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
    "use_bluesky": False,
    "GPT_summary": False,
    "broad_terms": ["e.g.", "multiple sclerosis", "MS", "autoimmunity"
                    ],
    "narrow_terms": ["e.g.", "T cells", "CD8", "TCR"
                     ],
    "DAYS_BACK": 7,
    "FILTER_SCOPE": "title",
    "prompt1": (
        "You are an expert summarizing recent scientific publications and posts. "
        "Focus on extracting key findings, trends, or notable advances relevant "
        "to the specified topic. Structure your summary as follows: write 2 sentences highlighting "
        "the main insights followed by a bullet point section for (if present) PubMed, bioRxiv, and Bluesky sources. "
        "If content is not relevant to the topic, omit it. Be succinct and professional."
    )
}
