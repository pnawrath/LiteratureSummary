from flask import Flask, render_template, request, jsonify
from config_file import BLUESKY_FEEDS, BIORXIV_SUBJECTS
import core
import markdown
import tracemalloc
from threading import Thread
import uuid
import re
import random

app = Flask(__name__)
tracemalloc.start()
task_store = {}  # task_id -> {status, result}


def run_summary_task(task_id, all_results, grouped_results, config):
    try:
        result = core.summarize_results(all_results, grouped_results, config).strip()
        task_store[task_id] = {"status": "done", "result": result}
    except Exception as e:
        task_store[task_id] = {"status": "error", "result": str(e)}


def extract_terms_from_gpt_output(text):
    lines = re.split(r"\n|,|;", text)
    terms = []
    for line in lines:
        cleaned = re.sub(r"^\s*\d+[\.\)]?\s*", "", line)
        if cleaned.strip():
            terms.append(cleaned.strip())
    return terms


@app.route('/status/<task_id>')
def task_status(task_id):
    return jsonify(task_store.get(task_id, {"status": "unknown", "result": ""}))


@app.route('/', methods=['GET', 'POST'])
def index():
    config = core.default_config.copy()
    summary = ""
    status_messages = []
    task_id = None

    subject = request.form.get("biorxiv_subject", "Immunology")
    config["biorxiv_subject"] = subject
    config["feed_uri"] = request.form.get("feed_uri", config["feed_uri"])
    html_help_text = markdown.markdown(config["HELP_TEXT"])

    if request.method == 'POST':
        interest_sentence = request.form.get("interest_sentence", "").strip()

        # Parse checkboxes here (common)
        config["GPT_summary"] = (request.form.get("GPT_summary") == "on")
        config["use_pubmed"] = (request.form.get("use_pubmed") == "on")
        config["use_biorxiv"] = (request.form.get("use_biorxiv") == "on")
        config["use_bluesky"] = (request.form.get("use_bluesky") == "on")

        if interest_sentence:
            config["interest_sentence"] = interest_sentence

            # Generate terms
            gpt_output = core.generate_terms_from_interest(interest_sentence, config)
            all_terms = extract_terms_from_gpt_output(gpt_output)
            random.shuffle(all_terms)

            # Split terms into two halves
            half = len(all_terms) // 2
            terms1, terms2 = all_terms[:half], all_terms[half:]

            # Use terms as broad and narrow directly (only one fetch)
            config["broad_terms"] = terms1
            config["narrow_terms"] = terms2

            all_results, grouped_results, counts, searched_sources = core.fetch_all_results(config)

        else:
            # Use existing form values
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

            all_results, grouped_results, counts, searched_sources = core.fetch_all_results(config)

        total_results_count = sum(counts.values())

        if total_results_count == 0:
            status_messages.append("No results found for your search terms in the selected sources.")

        # Launch summary task if enabled
        if config.get("GPT_summary") and all_results:
            if total_results_count < 60:
                task_id = str(uuid.uuid4())
                task_store[task_id] = {"status": "pending", "result": ""}
                thread = Thread(target=run_summary_task, args=(task_id, all_results, grouped_results, config))
                thread.start()
            else:
                summary = "GPT summary skipped because more than 60 results were found."

        current, peak = tracemalloc.get_traced_memory()
        print(f"Current memory usage: {current / 1024:.2f} KB")
        print(f"Peak memory usage: {peak / 1024:.2f} KB")

    else:
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
        task_id=task_id,
    )


if __name__ == "__main__":
    app.run(debug=True)