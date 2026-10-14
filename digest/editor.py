"""Claude steps, run through the Claude Code CLI on this machine (`claude -p`),
so they use the Claude Code login instead of a paid API key.

1. triage: score every new paper against interests.md
2. write_issue: pick and describe the week's papers from the shortlist

Links and titles are filled in by the renderer, never by the model, so the
issue cannot point at a paper that does not exist."""
import json
import re
import shutil
import subprocess
from typing import List, Literal

from pydantic import BaseModel, Field, ValidationError

Desk = Literal["tcr", "methods", "ebv_ms", "big_story", "off_topic", "other"]


class Score(BaseModel):
    id: int
    score: int = Field(description="0-10 relevance for this reader")
    desk: Desk


class Triage(BaseModel):
    scores: List[Score]


class Pick(BaseModel):
    id: int
    line: str = Field(description="One plain sentence, or empty")


class Issue(BaseModel):
    editorial: str = Field(description="One or two plain sentences on the week, or empty")
    front_page: List[Pick] = Field(description="The 3-5 papers this reader must not miss")
    preprints: List[Pick]
    ebv_ms: List[Pick]
    tcr: List[Pick]
    methods: List[Pick]
    social_media: List[str] = Field(description="2-5 one-line bullets on what the field is discussing")
    notable: List[Pick] = Field(description="Worth a glance; line may be empty")
    off_topic: List[Pick] = Field(description="1-3 general-interest pieces; line may be empty")


TRIAGE_RULES = """You triage new literature for one scientist. Their interests sheet is below.

Score each item 0-10 for how much THIS reader would want to see it:
- 9-10: directly on a core interest, or a major result in their field
- 6-8: clearly relevant, solid work
- 3-5: tangential
- 0-2: off-topic or low value

Desks: tcr = TCR specificity, cross-reactivity, structure, repertoire; methods = new
experimental or computational T cell methods; ebv_ms = EBV and multiple sclerosis;
big_story = a major immunology result of broad interest regardless of topic;
off_topic = a news feature, essay or science-writing piece on a general topic (not a
research paper), scored by how interesting it is to a curious scientist; other = anything else.

Be strict: most items each week should score below 5. Apply the reader's exclusions
literally. Return a score for every id given."""

STYLE = """Write the way a scientist talks to a colleague over coffee. Plain words, short sentences.
- Say what was done and what came out, using the concrete nouns and numbers in the abstract.
  If the abstract gives no specifics, keep the line short or leave it empty.
- Add why it matters to this reader only when that is not obvious, in a few words.
- One sentence per entry, ideally under 20 words.
- No intensifiers or sales words: genuinely, truly, really, remarkably, particularly,
  importantly, notably, new, newly, novel, first-ever, standout, must-read, groundbreaking, crucial, pivotal, key, exciting,
  promising, robust, comprehensive, landscape, underscores, highlights, sheds light,
  paves the way, delve, unveils.
- No contrast framing: never "not just X", "not another X", "more than X", "X, not Y",
  "rather than", "beyond X". State what the thing is, not what it is not.
- No colons for effect, no em dashes, no exclamation marks, no rhetorical questions.
- Do not start with "This paper", "This study", "The authors", "Researchers"."""

EDITOR_RULES = f"""You edit a short weekly literature digest for one scientist (interests below).
From the shortlisted items and social media posts, choose and describe this week's items.

{STYLE}

Sections (an id appears in at most one section; empty sections are fine):
- front_page: the 3-5 most important items for this reader, preprints included.
- preprints: every other selected preprint, whatever its topic.
- ebv_ms, tcr, methods: selected journal papers by topic.
- social_media: what people in the field are saying this week. Name the paper or topic
  and what the reaction was. Never mention the shortlist, the posts as a dataset, or counts., praising or criticising. Skip if nothing stands out.
- notable: other journal papers worth a glance, including big immunology stories.
- off_topic: 1-3 general-interest pieces (news features, essays) if any are good.
- editorial: one or two sentences about the week, or empty if there is nothing concrete to say.
Use only ids from the shortlist. Never invent papers, numbers or authors."""

BANNED = re.compile(
    r"genuine|\btruly|\breally|remarkabl|particularly|importantly|notabl|\bnovel|first.ever|"
    r"\bnew(ly|est)?\b|standout|must.read|long.awaited|\bfinally\b|at last|"
    r"groundbreaking|crucial|pivotal|exciting|promising|\brobust|comprehensive|landscape|"
    r"underscor|highlight|sheds? light|paves? the way|delve|unveil|game.changer|holistic|"
    r"\bnot just\b|\bnot another\b|\bmore than (just|a)\b|\brather than\b|\bbeyond\b|"
    r"^this (paper|study)|^the authors|^researchers|\u2014|!", re.IGNORECASE | re.MULTILINE)


def _schema(model):
    """Pydantic schema with $refs inlined, which the CLI accepts most reliably."""
    schema = model.model_json_schema()
    defs = schema.pop("$defs", {})

    def inline(node):
        if isinstance(node, dict):
            if "$ref" in node:
                return inline(defs[node["$ref"].split("/")[-1]])
            out = {k: inline(v) for k, v in node.items() if k != "title"}
            if out.get("type") == "object":
                out["additionalProperties"] = False
            return out
        if isinstance(node, list):
            return [inline(v) for v in node]
        return node

    return inline(schema)


class Editor:
    def __init__(self, cfg, interests):
        self.cfg = cfg["claude"]
        self.interests = interests
        self.binary = self.cfg.get("binary") or shutil.which("claude") or "claude"

    def _call(self, rules, content, model_cls, effort):
        cmd = [self.binary, "-p",
               "--output-format", "json",
               "--json-schema", json.dumps(_schema(model_cls)),
               "--system-prompt", f"{rules}\n\n<interests>\n{self.interests}\n</interests>",
               "--tools", "",
               "--no-session-persistence",
               "--model", self.cfg.get("model", "opus"),
               "--effort", effort]
        for attempt in range(2):
            proc = subprocess.run(cmd, input=content, capture_output=True, text=True,
                                  timeout=self.cfg.get("timeout_s", 900))
            try:
                result = json.loads(proc.stdout)
            except json.JSONDecodeError:
                result = {"is_error": True, "result": (proc.stderr or proc.stdout)[-500:]}
            if not result.get("is_error") and result.get("structured_output") is not None:
                try:
                    return model_cls.model_validate(result["structured_output"])
                except ValidationError as e:
                    error = f"invalid output ({e.error_count()} errors)"
            else:
                error = str(result.get("result") or result.get("subtype") or "unknown error")
        raise RuntimeError(f"Claude Code call failed: {error}")

    def triage(self, papers):
        size = self.cfg.get("triage_batch_size", 60)
        scores = {}
        for i in range(0, len(papers), size):
            batch = papers[i:i + size]
            lines = []
            for j, p in enumerate(batch, start=i):
                tag = " [general journal]" if p.get("general") else ""
                lines.append(f"[{j}] {p['title']} | {p['journal']}{tag} | {p['kind']}\n"
                             f"{p['abstract'][:1500] or '(no abstract)'}")
            out = self._call(TRIAGE_RULES, "\n\n".join(lines), Triage,
                             self.cfg.get("triage_effort", "low"))
            for s in out.scores:
                if i <= s.id < i + len(batch):
                    scores[s.id] = s
        return scores

    def write_issue(self, shortlist, posts):
        lines = []
        for p in shortlist:
            buzz = f" | discussed in {p['_buzz']} posts" if p.get("_buzz") else ""
            lines.append(f"[{p['_id']}] {p['title']} | {p['journal']} | {p['kind']} | "
                         f"desk guess: {p['_desk']}{buzz}\n{p['abstract'][:2500] or '(no abstract)'}")
        top_posts = sorted(posts, key=lambda p: p["likes"] + 2 * p["reposts"], reverse=True)[:150]
        post_lines = [f"- @{p['author']} ({p['likes']} likes): {p['text']}"
                      + (f" [link: {p['link_title']}]" if p["link_title"] else "") for p in top_posts]
        content = ("<shortlist>\n" + "\n\n".join(lines) + "\n</shortlist>\n\n"
                   "<bluesky>\n" + ("\n".join(post_lines) or "(no posts)") + "\n</bluesky>")
        effort = self.cfg.get("editor_effort", "high")
        issue = self._call(EDITOR_RULES, content, Issue, effort)

        # One rewrite pass if banned phrasing slipped through.
        offenders = sorted({m.group(0).lower() for t in _texts(issue) for m in BANNED.finditer(t)})
        if offenders:
            fix = (f"{content}\n\n<draft>\n{issue.model_dump_json()}\n</draft>\n\n"
                   f"Rewrite the draft without these words or patterns: {', '.join(offenders)}. "
                   "Keep the same ids and structure.")
            issue = self._call(EDITOR_RULES, fix, Issue, effort)
        return issue


def _texts(issue):
    yield issue.editorial
    for section in (issue.front_page, issue.preprints, issue.ebv_ms, issue.tcr,
                    issue.methods, issue.notable, issue.off_topic):
        for pick in section:
            yield pick.line
    yield from issue.social_media
