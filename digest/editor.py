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

Desk = Literal["tcr", "methods", "ebv_ms", "big_story", "other"]


class Score(BaseModel):
    id: int
    score: int = Field(description="0-10 relevance for this reader")
    desk: Desk


class Triage(BaseModel):
    scores: List[Score]


class Pick(BaseModel):
    id: int
    line: str = Field(description="One plain sentence: what was found, and why this reader cares")


class Issue(BaseModel):
    editorial: str = Field(description="Two or three plain sentences on what moved this week")
    front_page: List[Pick] = Field(description="The 3-5 papers this reader must not miss")
    big_stories: List[Pick]
    tcr: List[Pick]
    methods: List[Pick]
    ebv_ms: List[Pick]
    bluesky: List[str] = Field(description="2-5 one-line bullets on what the field is discussing")
    also_notable: List[int] = Field(description="ids worth a glance, no text")


TRIAGE_RULES = """You triage new immunology literature for one scientist. Their interests sheet is below.

Score each paper 0-10 for how much THIS reader would want to see it:
- 9-10: directly on a core interest, or a field-changing immunology story
- 6-8: clearly relevant, solid work
- 3-5: tangential
- 0-2: off-topic or low value

Desks: tcr = TCR specificity, cross-reactivity, structure, repertoire; methods = new
experimental or computational T cell methods; ebv_ms = EBV and multiple sclerosis
(especially CD8 T cells); big_story = a major immunology result of broad interest
regardless of topic (usually top journals); other = anything else.

Be strict: most papers each week should score below 5. Apply the reader's exclusions
literally. Return a score for every id given."""

STYLE = """Write like a sharp colleague summarising the week at lab meeting, not like a press release.
- Plain declarative sentences. Lead with the finding: what was shown, in what system, how.
- Concrete nouns and numbers from the abstract. If the abstract has no specifics, say less.
- Then, briefly, why it matters to this reader, only if that is not obvious.
- One sentence per paper, ideally under 25 words.
- Never use: groundbreaking, novel, crucial, pivotal, landscape, underscores, highlights,
  sheds light, paves the way, delve, intricate, robust, comprehensive, notably, unveils,
  game-changer, exciting, promising, holistic, "this study", "the authors".
- No em dashes, no exclamation marks, no rhetorical questions, no hedging filler."""

EDITOR_RULES = f"""You edit a short weekly literature digest for one immunologist (interests below).
From the shortlisted papers and Bluesky posts, choose and describe this week's papers.

{STYLE}

- front_page: the 3-5 most important papers for this reader. Do not repeat them elsewhere.
- Each desk: only papers that earn a place; an empty desk is fine.
- bluesky: what people are actually discussing, praising or criticising; skip if nothing stands out.
- editorial: two or three sentences, same style, no throat-clearing.
- also_notable: remaining shortlisted ids worth a glance.
- Use only ids from the shortlist. Never invent papers, numbers or authors."""

BANNED = re.compile(
    r"groundbreaking|\bnovel\b|crucial|pivotal|landscape|underscor|highlight|sheds? light|"
    r"paves? the way|delve|intricate|\brobust|comprehensive|notably|unveil|game.changer|"
    r"exciting|promising|holistic|this study|the authors|—", re.IGNORECASE)


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
                error = str(result.get("result") or result.get("subtype") or "unknown error")[:300]
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
        lines = [f"[{p['_id']}] {p['title']} | {p['journal']} | {p['kind']} | desk guess: {p['_desk']}\n"
                 f"{p['abstract'][:2500] or '(no abstract)'}" for p in shortlist]
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
    for section in (issue.front_page, issue.big_stories, issue.tcr, issue.methods, issue.ebv_ms):
        for pick in section:
            yield pick.line
    yield from issue.bluesky
