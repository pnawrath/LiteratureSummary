from flask import Flask, render_template, request
import core

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def index():
    config = core.default_config.copy()  # start with defaults
    summary = ""
    status_messages = []  # Keep for any status or errors

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

        # Optional: add some status messages if needed
        # For example, if no results found anywhere:
        if all(len(grouped_results[src]) == 0 for src in searched_sources):
            status_messages.append("No results found for your search terms in the selected sources.")

        # Generate GPT summary if requested and if we have results
        if config.get("GPT_summary") and all_results:
            summary = core.summarize_results(all_results, config).strip()

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
        searched_sources=searched_sources
    )


if __name__ == "__main__":
    app.run(debug=True)
