#!/usr/bin/env bash
# SY Gate 1c driver —— Tier A SNP 去冗 + sanity 阳性对照
#
# DO NOT export $HOME/envs/snakemake/bin onto PATH here. Doing so makes the
# snakemake env's python3.14 (no numpy) win over the per-rule conda env that
# snakemake activates, and every script: rule that imports numpy dies with
# ModuleNotFoundError while the env itself is perfectly fine.
# snakemake is always invoked by absolute path instead.
set -uo pipefail
TASK="${1:-sy_gate1c}"
TASKDIR="$HOME/work/SY"
R="$TASKDIR/run/$TASK"
PDIR="$TASKDIR/pipeline/plague-decay"
SP="$HOME/envs/snakemake/bin/snakemake"
CP="$HOME/envs/snakemake-envs"
mkdir -p "$TASKDIR/logs" "$R/steps"
date +%s > "$R/start"
trap 'echo "$? $(date +%s)" > "$R/exit"' EXIT
export PYTHONUNBUFFERED=1
export CONDA_PKGS_DIRS="$HOME/.conda/pkgs"
export CONDA_ENVS_PATH="$CP"

# a python that actually has numpy, for the out-of-DAG sanity check
PYPY=""
for p in "$CP"/*/bin/python; do
  if "$p" -c "import numpy" >/dev/null 2>&1; then PYPY="$p"; break; fi
done

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

echo "=== gate1c start $(date) ; PYPY=$PYPY ==="
cd "$PDIR" || { echo "CDFAIL $PDIR"; exit 2; }

step gate1c 14400 bash -c "$SP --use-conda --conda-prefix $CP --cores 32 --rerun-incomplete --latency-wait 60 --keep-going gate1c"

step sanity 300 bash -c "
  if [ -n '$PYPY' ]; then
    '$PYPY' scripts/gate1c_sanity.py \
      results/dedup/snp_distance_tierA.npy \
      results/dedup/genotypes_tierA.npz \
      results/dedup/tierA_clusters.tsv \
      results/dedup/tierA_keep.tsv; echo SANITY_EXIT=\$?
  else
    echo 'no numpy python found'; echo SANITY_EXIT=99
  fi
"

step report 180 bash -c '
  echo "=== Tier A keep ==="
  test -f results/dedup/tierA_keep.tsv && echo "  n=$(( $(wc -l < results/dedup/tierA_keep.tsv) - 1 ))"
  for L in 02_genotype_matrix_tierA 02_snp_distance_tierA 02_tierA_dedup; do
    echo "--- $L ---"; cat logs/$L.log 2>/dev/null | head -14; done
'
echo "=== gate1c done $(date) ==="
