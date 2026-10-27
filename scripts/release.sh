#!/usr/bin/env bash
# Local release (no GitHub Actions, so no 'workflow' token scope needed).
# Usage: scripts/release.sh 1.2.1      -> tests (also without matplotlib/openpyxl), builds dist/b2s.pyz,
#                                         tags v1.2.1, pushes, creates the GitHub release with b2s.pyz attached.
set -euo pipefail
V=${1:?version, e.g. 1.2.1}
cd "$(dirname "$0")/.."
PY=${PY:-python3}
$PY -m pytest -q
$PY scripts/build_pyz.py
$PY -I dist/b2s.pyz --version | grep -q "$V" || { echo "b2s.pyz version != $V (update __init__.py)"; exit 1; }
git tag -a "v$V" -m "busco2summary $V"
for i in 1 2 3 4 5; do git push origin main "v$V" && break; sleep 5; done   # VPN proxy sometimes times out
for i in 1 2 3 4 5; do
  gh release create "v$V" dist/b2s.pyz --title "busco2summary $V" --generate-notes && break
  sleep 5; gh release view "v$V" >/dev/null 2>&1 && break
done
gh release view "v$V" --json url -q .url
