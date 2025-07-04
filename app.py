from flask import Flask, render_template, request
import core

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def index():
    config = core.default_config.copy()  # start with defaults
    results = []
    summary = ""
    status_messages = []  # Ensure always defined

    if request.method == 'POST':
        for key in config:
            val = request.form.get(key)
            if val is not None:
                if key in ["use_pubmed", "use_biorxiv", "use_bluesky", "GPT_summary"]:
                    config[key] = (val == "on")
                elif key == "DAYS_BACK":
                    config[key] = int(val)
                elif key in ["broad_terms", "narrow_terms"]:
                    config[key] = [t.strip() for t in val.split(",")]
                else:
                    config[key] = val
            else:
                if key in ["use_pubmed", "use_biorxiv", "use_bluesky", "GPT_summary"]:
                    config[key] = False

        if config["use_pubmed"]:
            pubmed_results, pubmed_status = core.fetch_pubmed_results(config)
            results.extend(pubmed_results)
            status_messages.append(pubmed_status)
        if config["use_biorxiv"]:
            results.extend(core.fetch_biorxiv_results(config))
        if config["use_bluesky"]:
            results.extend(core.fetch_custom_bluesky_feed(config))
        if config["GPT_summary"] and results:
            summary = core.summarize_results(results, config).strip()

    # Always define grouped_results
    grouped_results = {
        "PubMed": [],
        "bioRxiv": [],
        "Bluesky": []
    }
    for r in results:
        src = r.get("source", "")
        if src in grouped_results:
            grouped_results[src].append(r)

    return render_template(
        'index.html',
        config=config,
        results=results,
        summary=summary,
        grouped_results=grouped_results,
        status_messages=status_messages
    )


if __name__ == "__main__":
    app.run(debug=True)
