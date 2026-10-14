"""Turns the editor's Issue into an Obsidian note."""
from collections import Counter


def _entry(p, line=""):
    """Every entry, in every section: [Title](url) · *Journal* · Author et al. · Text"""
    title = p["title"].rstrip(".").replace("[", "(").replace("]", ")")
    parts = [f"[{title}]({p['url']})" if p["url"] else title, f"*{p['journal']}*"]
    if p["authors"]:
        parts.append(p["authors"])
    if line and line.strip():
        parts.append(line.strip())
    return "- " + " · ".join(parts)


SECTIONS = (("Front Page", "front_page"), ("Preprints", "preprints"), ("EBV & MS", "ebv_ms"),
            ("TCR specificity & structure", "tcr"), ("Methods", "methods"))
TAIL = (("Notable", "notable"), ("Off Topic", "off_topic"))


def render(issue, by_id, week, start, end, stats, problems):
    num = week.split("-W")[1].lstrip("0")
    out = [
        "---",
        "type: literature-digest",
        f"week: {week}",
        f"period: {start.isoformat()} to {end.isoformat()}",
        f"screened: {stats['screened']}",
        "tags: [literature-digest]",
        "---",
        f"# Week {num} · {start.day} {start:%b} to {end.day} {end:%b %Y}",
        "",
    ]
    if issue.editorial.strip():
        out += [f"> {issue.editorial.strip()}", ""]

    used = set()

    def section(heading, picks):
        rows = []
        for pk in picks:
            if pk.id in by_id and pk.id not in used:
                used.add(pk.id)
                rows.append(_entry(by_id[pk.id], pk.line))
        return [f"## {heading}", *rows, ""] if rows else []

    for heading, field in SECTIONS:
        out += section(heading, getattr(issue, field))
    if issue.social_media:
        out += ["## Social Media", *[f"- {b}" for b in issue.social_media], ""]
    for heading, field in TAIL:
        out += section(heading, getattr(issue, field))

    per_source = Counter(stats["sources"])
    src = ", ".join(f"{k} {v}" for k, v in per_source.most_common())
    out += ["---",
            f"*Screened {stats['screened']} new papers ({src}) and {stats['posts']} Bluesky posts; "
            f"{stats['skipped_seen']} already seen in earlier issues were skipped.*"]
    if problems:
        out += ["", "> [!warning] Source problems this week",
                *[f"> - {p}" for p in problems]]
    return "\n".join(out) + "\n"


def update_index(index_path, week, rel_link, headline):
    header = "# Literature Digest\n\nNewest issue first.\n\n"
    lines = []
    if index_path.exists():
        lines = [l for l in index_path.read_text().splitlines() if l.startswith("- [[")]
    lines = [l for l in lines if f"|{week}]]" not in l]
    lines.insert(0, f"- [[{rel_link}|{week}]] · {headline}")
    index_path.write_text(header + "\n".join(lines) + "\n")
