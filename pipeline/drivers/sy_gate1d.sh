#!/usr/bin/env bash
# SY Gate 1d —— Tier B: 3,441 downloads + snippy + backbone dedup (10-14 h)
# 注意：不要把 ~/envs/snakemake/bin 放进 PATH（会遮蔽 per-rule env 的 python）
set -uo pipefail
TASK="${1:-sy_gate1d}"
TASKDIR="$HOME/work/SY"; R="$TASKDIR/run/$TASK"
PDIR="$TASKDIR/pipeline/plague-decay"
SP="$HOME/envs/snakemake/bin/snakemake"; CP="$HOME/envs/snakemake-envs"
mkdir -p "$TASKDIR/logs" "$R/steps"
date +%s > "$R/start"
trap 'echo "$? $(date +%s)" > "$R/exit"' EXIT
export PYTHONUNBUFFERED=1 CONDA_PKGS_DIRS="$HOME/.conda/pkgs" CONDA_ENVS_PATH="$CP"

step() {
  local n="$1" t="$2"; shift 2
  [[ -f "$R/steps/$n.done" ]] && { echo "[skip] $n"; return 0; }
  echo "[step:$n] $(date +%H:%M:%S) start"
  date +%s > "$R/steps/$n.start"
  timeout -k 60 "$t" "$@"; local rc=$?
  [[ $rc -ne 0 ]] && { echo "[step:$n] FAILED rc=$rc"; return $rc; }
  date +%s > "$R/steps/$n.done"; echo "[step:$n] $(date +%H:%M:%S) done"
}

echo "=== gate1d start $(date) ==="
cd "$PDIR" || exit 2
step gate1d 345600 bash -c "$SP --use-conda --conda-prefix $CP --cores 32 --rerun-incomplete --latency-wait 60 --keep-going gate1d"
step report 300 bash -c '
  echo "=== Tier B backbone ==="
  test -f results/dedup/tierB_backbone.tsv && echo "  backbone n=$(( $(wc -l < results/dedup/tierB_backbone.tsv) - 1 ))"
  cat logs/02_tierB_backbone.log 2>/dev/null | head -14
  echo "=== blocked pairs (校准 min_joint_called_frac 用) ==="
  cat logs/02_snp_distance_tierB.log 2>/dev/null | grep -iE "blocked|joint|percentile" | head -12
'
echo "=== gate1d done $(date) ==="
