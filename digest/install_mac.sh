#!/bin/zsh
# Usage: digest/install_mac.sh "/path/to/your/Obsidian vault"
# Sets up a Python environment, the Literature Digest folder and the weekday schedule.
set -e
VAULT="${1:?Pass the path to your Obsidian vault}"
DIGEST="$(cd "$(dirname "$0")" && pwd)"
FOLDER="$VAULT/Literature Digest"
LOG="$HOME/Library/Logs/literature-digest.log"
PLIST="$HOME/Library/LaunchAgents/com.literaturedigest.weekly.plist"

command -v claude >/dev/null || { echo "Claude Code (claude) not found on PATH. Install it and log in first."; exit 1; }

python3 -m venv "$DIGEST/.venv"
"$DIGEST/.venv/bin/pip" install -q -r "$DIGEST/requirements.txt"

mkdir -p "$FOLDER/Issues"
[ -f "$FOLDER/interests.md" ] || cp "$DIGEST/interests.example.md" "$FOLDER/interests.md"
[ -f "$DIGEST/.env" ] || echo "ENTREZ_EMAIL=" > "$DIGEST/.env"

sed "s|__VAULT__|$VAULT|" "$DIGEST/launchd/run_digest.sh" > "$DIGEST/launchd/run_digest.local.sh"
chmod +x "$DIGEST/launchd/run_digest.local.sh"
mkdir -p "$HOME/Library/LaunchAgents"
sed -e "s|__SCRIPT__|$DIGEST/launchd/run_digest.local.sh|" -e "s|__LOG__|$LOG|" \
    "$DIGEST/launchd/com.literaturedigest.weekly.plist" > "$PLIST"
launchctl bootout "gui/$(id -u)" "$PLIST" 2>/dev/null || true
launchctl bootstrap "gui/$(id -u)" "$PLIST"

echo "Installed. Issues go to: $FOLDER/Issues"
echo "Edit $FOLDER/interests.md and add your email to $DIGEST/.env (NCBI asks for one)."
echo "Run now:   launchctl kickstart gui/$(id -u)/com.literaturedigest.weekly"
echo "Log:       $LOG"
