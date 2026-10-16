#!/bin/zsh
# Called by launchd on weekdays at 10:30. Exits quickly if this week's issue exists.
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
export PATH="$HOME/.local/bin:$HOME/.claude/local:/opt/homebrew/bin:/usr/local/bin:$PATH"
cd "$REPO" || exit 1
echo "=== $(date) ==="
exec "$REPO/digest/.venv/bin/python" -m digest.run --vault "__VAULT__"
