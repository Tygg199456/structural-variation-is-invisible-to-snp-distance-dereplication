#!/usr/bin/env bash
# SY: Tier A dedup under three profile-constraint modes + mixed-cluster
# classification (external review 2026-10-11, execution steps 2-3-4).
set -uo pipefail

TASK="${1:-sy_dedup_modes}"; TASKDIR="$HOME/work/SY"; R="$TASKDIR/run/$TASK"
PDIR="$TASKDIR/pipeline/plague-decay"; CP="$HOME/envs/snakemake-envs"

mkdir -p "$R/steps" "$TASKDIR/logs"
date +%s > "$R/start"
trap 'echo "$? $(date +%s)" > "$R/exit"' EXIT
export PYTHONUNBUFFERED=1

step() {
  local n="$1" t="$2"; shift 2
  [[ -f "$R/steps/$n.done" ]] && { echo "[skip] $n"; return 0; }
  echo "[step:$n] $(date +%H:%M:%S) start"
  date +%s > "$R/steps/$n.start"
  timeout -k 60 "$t" "$@"; local rc=$?
  [[ $rc -ne 0 ]] && { echo "[step:$n] FAILED rc=$rc"; return $rc; }
  date +%s > "$R/steps/$n.done"
  echo "[step:$n] $(date +%H:%M:%S) done"
}

cd "$PDIR" || { echo "CDFAIL $PDIR"; exit 2; }

PYPY=""
for p in "$CP"/*/bin/python; do
  "$p" -c "import numpy,scipy" 2>/dev/null && { PYPY="$p"; break; }
done
[[ -n "$PYPY" ]] || { echo "FATAL: no python with numpy+scipy under $CP"; exit 2; }
echo "PYPY=$PYPY"
echo "=== sy_dedup_modes start $(date) ==="

step modes 3600 bash -c "set -o pipefail; \"$PYPY\" scripts/dedup_modes.py 2>&1 | tee results/qc/dedup_modes.txt" || exit 1

# The production keep/clusters must be unchanged (mode default = full). Assert
# the full-mode run reproduces the shipped 206.
step assert 120 bash -c '
  n=$(wc -l < results/dedup/tierA_keep.full.tsv)
  echo "full-mode keep rows (incl header) = $n"
  [ "$n" -eq 207 ] || { echo "FATAL: full mode did not reproduce 206"; exit 1; }
  echo "OK: full mode reproduces the shipped 206"
' || exit 1

echo "=== sy_dedup_modes done $(date) ==="
