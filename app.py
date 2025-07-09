from flask import Flask, render_template, request
from config_file import BLUESKY_FEEDS, BIORXIV_SUBJECTS, HELP_TEXT
import core
import markdown

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def index():
    config = core.default_config.copy()  # start with defaults
    summary = ""
    status_messages = []  # Keep for any status or errors
    subject = request.form.get("biorxiv_subject", "Immunology")
    config["biorxiv_subject"] = subject
    config["feed_uri"] = request.form.get("feed_uri", config["feed_uri"])
    html_help_text = markdown.markdown(config["HELP_TEXT"])

    if request.method == 'POST':
        # Update config from form data
        for key in config:
            val = request.form.get(key)
            if val is not None:
                if key in ["use_pubmed", "use_biorxiv", "use_bluesky", "GPT_summary"]:
                    config[key] = (val == "on")
                elif key == "DAYS_BACK":
                    try:
                        config[key] = int(val)
                    except ValueError:
                        config[key] = core.default_config.get(key, 30)
                elif key in ["broad_terms", "narrow_terms"]:
                    config[key] = [t.strip() for t in val.split(",") if t.strip()]
                else:
                    config[key] = val
            else:
                if key in ["use_pubmed", "use_biorxiv", "use_bluesky", "GPT_summary"]:
                    config[key] = False

        # Fetch results with unified function
        all_results, grouped_results, counts, searched_sources = core.fetch_all_results(config)
        total_results_count = sum(len(grouped_results[src]) for src in searched_sources)

        # Optional: add some status messages if needed
        # For example, if no results found anywhere:
        if all(len(grouped_results[src]) == 0 for src in searched_sources):
            status_messages.append("No results found for your search terms in the selected sources.")

        # Generate GPT summary if requested and if we have results
        summary = ""
        if config.get("GPT_summary") and all_results and total_results_count < 100:
            summary = core.summarize_results(all_results, config).strip()
        elif config.get("GPT_summary") and total_results_count >= 100:
            summary = "GPT summary skipped because more than 100 results were found."

    else:
        # GET request — empty defaults
        all_results = []
        grouped_results = {"PubMed": [], "bioRxiv": [], "Bluesky": []}
        counts = {"PubMed": 0, "bioRxiv": 0, "Bluesky": 0}
        searched_sources = []

    return render_template(
        'index.html',
        config=config,
        results=all_results,
        summary=summary,
        grouped_results=grouped_results,
        status_messages=status_messages,
        counts=counts,
        searched_sources=searched_sources,
        biorxiv_subjects=BIORXIV_SUBJECTS,
        bluesky_feeds=BLUESKY_FEEDS,
        help_text=html_help_text,
    )

if __name__ == "__main__":
    app.run(debug=True)
