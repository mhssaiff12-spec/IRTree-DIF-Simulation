#!/usr/bin/env bash
# Safe one-time creation of the results branch from the included backup.
# Does not delete/overwrite the source checkout.
set -euo pipefail
ROOT=$(git rev-parse --show-toplevel)
cd "$ROOT"
if git ls-remote --exit-code --heads origin results >/dev/null 2>&1; then
  echo "results branch already exists; not replacing it." >&2
  exit 1
fi
[ -d results_core8_backup ] || { echo "Missing results_core8_backup" >&2; exit 1; }
TMP=$(mktemp -d)
trap 'git worktree remove --force "$TMP/wt" >/dev/null 2>&1 || true; rm -rf "$TMP"' EXIT
git worktree add --detach "$TMP/wt" HEAD
cd "$TMP/wt"
git checkout --orphan results
git rm -r --cached . >/dev/null 2>&1 || true
mkdir -p sim/results_core8
cp -a "$ROOT/results_core8_backup/." sim/results_core8/
git add -f sim/results_core8
git -c user.name=irtree-bot -c user.email=irtree-bot@users.noreply.github.com commit -m 'Seed existing Claude replications'
git push origin HEAD:results
printf 'Results branch seeded safely.\n'
