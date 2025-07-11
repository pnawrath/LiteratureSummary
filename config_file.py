import os
from dotenv import load_dotenv

load_dotenv()
help_text_path = os.path.join("static", "info.md")
with open(help_text_path, "r", encoding="utf-8") as file:
    HELP_TEXT = file.read()

default_config = {
    "Entrez_email": os.getenv("ENTREZ_EMAIL", ""),
    "biorxiv_subject": "Immunology",
    "handle_bsky": os.getenv("HANDLE_BSKY", ""),
    "bsky_app_password": os.getenv("BSKY_APP_PASSWORD", ""),
    "feed_uri": os.getenv("FEED_URI", "at://did:plc:rwkarouaeku2g6qkvqkirwa5/app.bsky.feed.generator/ImmunoSky"),
    "openai_api_key": os.getenv("OPENAI_API_KEY", ""),
    "use_pubmed": True,
    "use_biorxiv": True,
    "use_bluesky": True,
    "GPT_summary": True,
    "broad_terms": [],
    "narrow_terms": [],
    "DAYS_BACK": 14,
    "interest_sentence": "",
    "FILTER_SCOPE": "title",
    "HELP_TEXT": HELP_TEXT,
    "prompt1": (
        "Focus on extracting key findings, trends, or notable advances. Structure your summary as follows: Write 2 "
        "brief sentences highlighting the main insights followed by a bullet-point section for (if present) PubMed, "
        "bioRxiv, and Bluesky. If content is not relevant to the topic, omit it. Be succinct and professional. "
    )
}

BLUESKY_FEEDS = {
    "T cells in MS": "at://did:plc:g3xk2d3jg6dma5glgwilyiyp/app.bsky.feed.generator/aaaknfoyzxnfa",
    "ImmunoSky": "at://did:plc:rwkarouaeku2g6qkvqkirwa5/app.bsky.feed.generator/ImmunoSky",
    "MLSky": "at://did:plc:rwkarouaeku2g6qkvqkirwa5/app.bsky.feed.generator/MLSky",
    "GeneSky": "at://did:plc:rwkarouaeku2g6qkvqkirwa5/app.bsky.feed.generator/GeneSky"
}

BIORXIV_SUBJECTS = [
    "Animal Behavior and Cognition", "Biochemistry", "Bioengineering", "Bioinformatics", "Biophysics",
    "Cancer Biology", "Cell Biology", "Clinical Trials", "Developmental Biology", "Ecology", "Epidemiology",
    "Evolutionary Biology", "Genetics", "Genomics", "Immunology", "Microbiology", "Molecular Biology",
    "Neuroscience", "Paleontology", "Pathology", "Pharmacology and Toxicology", "Physiology",
    "Plant Biology", "Scientific Communication and Education", "Synthetic Biology", "Systems Biology", "Zoology"
]