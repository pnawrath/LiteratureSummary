"""Turns the editor's Issue into an Obsidian note."""
from collections import Counter


def _ref(p):
    title = p["title"].rstrip(".").replace("[", "(").replace("]", ")")
    link = f"[{title}]({p['url']})" if p["url"] else title
    meta = " · ".join(x for x in (f"*{p['journal']}*", p["authors"]) if x)
    if p["kind"] == "preprint":
        meta += " · preprint"
    return f"{link} · {meta}"


def _section(heading, picks, by_id):
    rows = [f"- {_ref(by_id[pk.id])}  \n  {pk.line}" for pk in picks if pk.id in by_id]
    return [f"## {heading}", *rows, ""] if rows else []


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
        f"> {issue.editorial.strip()}",
        "",
    ]

    used = set()
    if issue.front_page:
        out.append("## Front page")
        for n, pk in enumerate(issue.front_page, 1):
            if pk.id in by_id:
                used.add(pk.id)
                out.append(f"{n}. **{_ref(by_id[pk.id])}**  \n   {pk.line}")
        out.append("")

    for heading, picks in (("Big immunology stories", issue.big_stories),
                           ("TCR specificity & structure", issue.tcr),
                           ("Methods & computation", issue.methods),
                           ("EBV & MS", issue.ebv_ms)):
        picks = [pk for pk in picks if pk.id not in used]
        used.update(pk.id for pk in picks)
        out += _section(heading, picks, by_id)

    if issue.bluesky:
        out += ["## On Bluesky", *[f"- {b}" for b in issue.bluesky], ""]

    notable = [i for i in issue.also_notable if i in by_id and i not in used]
    if notable:
        out += ["## Also notable", *[f"- {_ref(by_id[i])}" for i in notable], ""]

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
