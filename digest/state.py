"""Memory between issues: which papers were already considered, so nothing
shows up twice, and when the last issue ran."""
import json
import re
from datetime import date, timedelta
from pathlib import Path

KEEP_DAYS = 400


def title_key(title):
    return re.sub(r"[^a-z0-9]", "", title.lower())[:120]


def keys(item):
    out = [f"t:{title_key(item['title'])}"]
    if item.get("doi"):
        out.append(f"d:{item['doi']}")
    if item.get("pmid"):
        out.append(f"p:{item['pmid']}")
    return out


class State:
    def __init__(self, path: Path):
        self.path = path
        data = json.loads(path.read_text()) if path.exists() else {}
        self.seen = data.get("seen", {})          # key -> ISO date first seen
        self.last_run = data.get("last_run")      # ISO date of the last issue

    def is_seen(self, item):
        return any(k in self.seen for k in keys(item))

    def forget_since(self, day):
        """Undo everything recorded on or after `day`, so a rebuild of the
        current week sees the same papers again."""
        cutoff = day.isoformat()
        self.seen = {k: v for k, v in self.seen.items() if v < cutoff}
        if self.last_run and self.last_run >= cutoff:
            self.last_run = None

    def mark(self, items, today):
        for it in items:
            for k in keys(it):
                self.seen.setdefault(k, today.isoformat())

    def save(self, today):
        cutoff = (today - timedelta(days=KEEP_DAYS)).isoformat()
        self.seen = {k: v for k, v in self.seen.items() if v >= cutoff}
        self.last_run = today.isoformat()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps({"last_run": self.last_run, "seen": self.seen}, indent=0))
        tmp.replace(self.path)


def merge_duplicates(items):
    """Collapse the same paper arriving from several sources (e.g. PubMed and a
    journal's feed). Keeps the richest record and the most specific journal."""
    merged, index = [], {}
    for it in items:
        hit = next((index[k] for k in keys(it) if k in index), None)
        if hit is None:
            index.update({k: len(merged) for k in keys(it)})
            merged.append(dict(it))
            continue
        cur = merged[hit]
        if len(it["abstract"]) > len(cur["abstract"]):
            cur["abstract"] = it["abstract"]
        if cur["source"] == "PubMed" and it["source"] != "PubMed":
            cur["journal"], cur["general"] = it["journal"], it["general"]
        for f in ("doi", "pmid", "authors", "url"):
            cur[f] = cur[f] or it[f]
        index.update({k: hit for k in keys(it)})
    return merged


def drop_seen(items, state):
    fresh = [it for it in items if not state.is_seen(it)]
    return fresh, len(items) - len(fresh)


def count_mentions(item, posts):
    """How many posts link to this paper. Matched on DOI, which is the one
    identifier that survives a publisher's redirects."""
    doi = item.get("doi")
    if not doi:
        return 0
    return sum(1 for p in posts if doi in (p.get("link") or "").lower()
               or doi in (p.get("text") or "").lower())


def iso_week(today: date):
    y, w, _ = today.isocalendar()
    return f"{y}-W{w:02d}"
