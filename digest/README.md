# Literature Digest

A weekly literature digest that lands in your Obsidian vault as a short newspaper.
It uses the same sources as the web app (PubMed, bioRxiv, Bluesky) and adds journal tables of contents.

- **No API bill.** The Claude steps run through Claude Code (`claude -p`) with your Claude login.
- **Private.** It reads only `<vault>/Literature Digest/interests.md` and writes only inside that folder.
- **No repeats.** `.state/seen.json` remembers every paper already considered.
- **Hands-off.** A macOS LaunchAgent runs on weekdays at 10:30 and writes one issue per ISO week.
  On Monday it normally writes the issue; if the Mac was off, the next weekday catches up.

## Layout in the vault

```
Literature Digest/
  interests.md            <- you edit this
  Index.md                <- list of issues, newest first
  Issues/2026/2026-W42.md <- one per week
  .state/seen.json        <- memory (hidden from Obsidian)
```

## How an issue is built

1. Fetch the window since the last issue: PubMed queries, the bioRxiv immunology category,
   journal feeds and Bluesky feeds (all in `config.yaml`).
2. Merge duplicates across sources and drop papers seen in earlier issues.
3. Claude scores every paper 0 to 10 against `interests.md` (low effort, batches of 60).
4. The top 30 go to an editor step with a strict plain-language style guide; if banned
   phrasing slips through, the draft is rewritten once.
5. Titles, journals and links come from the source data, never from the model.

## Install (macOS)

Prerequisite: Claude Code installed and logged in (`claude` works in Terminal).

```
git clone https://github.com/pnawrath/LiteratureSummary && cd LiteratureSummary
digest/install_mac.sh "/path/to/your/Obsidian vault"
```

Then add your email to `digest/.env` (`ENTREZ_EMAIL=...`, NCBI asks for one), edit
`interests.md`, and test:

```
digest/.venv/bin/python -m digest.run --vault "/path/to/vault" --sources-only   # feeds only
launchctl kickstart gui/$(id -u)/com.literaturedigest.weekly                    # full run now
tail -f ~/Library/Logs/literature-digest.log
```

`--force` rebuilds the current week's issue.
