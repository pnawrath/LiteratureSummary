"""Build this week's Literature Digest issue in an Obsidian vault.

    python -m digest.run --vault ~/path/to/vault            # normal run
    python -m digest.run --vault ... --force                # rebuild this week's issue
    python -m digest.run --vault ... --sources-only         # test the feeds, no Claude

Reads only <vault>/<folder>/interests.md and writes only inside <folder>.
Safe to run every day: it exits early once this week's issue exists."""
import argparse
import os
import subprocess
import sys
from collections import Counter
from datetime import date, timedelta
from pathlib import Path

import yaml
from dotenv import load_dotenv

from . import render, sources
from .editor import Editor
from .state import People, State, count_mentions, drop_seen, iso_week, merge_duplicates

HERE = Path(__file__).resolve().parent
load_dotenv(HERE / ".env")  # ENTREZ_EMAIL, optional BSKY_HANDLE / BSKY_APP_PASSWORD


def log(msg):
    print(msg, flush=True)


def notify(title, text):
    if sys.platform == "darwin":
        script = f'display notification "{text}" with title "{title}"'
        subprocess.run(["osascript", "-e", script], check=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vault", default=os.getenv("DIGEST_VAULT"), help="Obsidian vault folder")
    ap.add_argument("--config", default=str(HERE / "config.yaml"))
    ap.add_argument("--force", action="store_true", help="rebuild even if this week's issue exists")
    ap.add_argument("--sources-only", action="store_true", help="fetch and report counts, no Claude")
    args = ap.parse_args(argv)
    if not args.vault:
        ap.error("--vault (or DIGEST_VAULT) is required")

    cfg = yaml.safe_load(Path(args.config).read_text())
    root = Path(args.vault).expanduser() / cfg["folder"]
    today = date.today()
    week = iso_week(today)
    issue_rel = Path("Issues") / str(today.isocalendar()[0]) / f"{week}.md"
    issue_path = root / issue_rel

    if issue_path.exists() and not args.force and not args.sources_only:
        log(f"{week} already exists, nothing to do.")
        return 0

    interests_path = root / "interests.md"
    if not interests_path.exists():
        log(f"Missing {interests_path}. Copy digest/interests.example.md there first.")
        return 1
    interests = interests_path.read_text()

    state = State(root / ".state" / "seen.json")
    if args.force:
        state.forget_since(today - timedelta(days=today.weekday()))
    start, end = sources.window(cfg, state.last_run)
    log(f"Window {start} to {end}")

    bcfg = cfg.get("bluesky", {})
    handles = bcfg.get("follows_of") or []
    people = People(root / ".state" / "follows.json")
    early = []
    if handles and people.stale(handles, today):
        log(f"Refreshing the people list from the follows of {', '.join(handles)}...")
        found, problem = sources.fetch_follow_set(handles, bcfg.get("max_people", 2000))
        if found:
            people.update(found, handles, today)
        if problem:
            early.append(problem)
    if people.people:
        log(f"Reading {len(people.people)} people on Bluesky")

    papers, posts, problems = sources.fetch_everything(cfg, start, end, people.people)
    problems = early + problems
    counts = Counter(p["source"] for p in papers)
    log(f"Fetched {len(papers)} papers {dict(counts)} and {len(posts)} Bluesky posts")
    for p in problems:
        log(f"  problem: {p}")

    if posts and bcfg.get("dynamic_terms", True) and not args.sources_only:
        editor_for_terms = Editor(cfg, interests)
        try:
            terms = editor_for_terms.suggest_terms(posts)
        except Exception as e:
            terms, _ = [], log(f"  problem: term suggestion failed ({e})")
        if terms:
            log("Terms from this week's posts: " + ", ".join(terms))
            extra, problem = sources.fetch_pubmed_terms(cfg, terms, start, end)
            log(f"  {len(extra)} further papers from those terms")
            papers.extend(extra)
            if problem:
                problems.append(problem)

    papers = merge_duplicates(papers)
    fresh, skipped = drop_seen(papers, state)
    log(f"{len(fresh)} new after removing duplicates and {skipped} seen in earlier issues")
    if args.sources_only:
        return 0
    if not fresh and not posts:
        log("Nothing new; no issue written.")
        return 1

    editor = Editor(cfg, interests)
    log("Triage...")
    scores = editor.triage(fresh)
    ccfg = cfg["claude"]
    ranked = []
    for i, p in enumerate(fresh):
        s = scores.get(i)
        if s and s.score >= ccfg.get("min_score", 5):
            ranked.append(dict(p, _id=i, _score=s.score, _desk=s.desk))
    for p in ranked:
        p["_buzz"] = count_mentions(p, posts)
    ranked.sort(key=lambda p: (p["_score"], p["_buzz"]), reverse=True)
    shortlist = ranked[:ccfg.get("shortlist_size", 45)]
    log(f"Shortlisted {len(shortlist)} of {len(fresh)}")

    log("Writing issue...")
    issue = editor.write_issue(shortlist, posts)
    by_id = {p["_id"]: p for p in shortlist}
    stats = {"screened": len(fresh), "sources": [p["source"] for p in fresh],
             "posts": len(posts), "skipped_seen": skipped}
    text = render.render(issue, by_id, week, start, end, stats, problems)

    issue_path.parent.mkdir(parents=True, exist_ok=True)
    issue_path.write_text(text)
    headline = by_id[issue.front_page[0].id]["title"] if issue.front_page and issue.front_page[0].id in by_id else ""
    render.update_index(root / "Index.md", week, week, headline[:90])

    state.mark(fresh, today)
    state.save(today)
    log(f"Wrote {issue_path}")
    notify("Literature Digest", f"{week} is ready in Obsidian")
    return 0


if __name__ == "__main__":
    sys.exit(main())
