#!/usr/bin/env bash
set -uo pipefail
TASK="${1:-sy_snippy1}"; TASKDIR="$HOME/work/SY"; R="$TASKDIR/run/$TASK"
PDIR="$TASKDIR/pipeline/plague-decay"; SP="$HOME/envs/snakemake/bin/snakemake"; CP="$HOME/envs/snakemake-envs"
mkdir -p "$TASKDIR/logs" "$R/steps"; date +%s > "$R/start"
trap 'echo "$? $(date +%s)" > "$R/exit"' EXIT
export PYTHONUNBUFFERED=1 CONDA_PKGS_DIRS="$HOME/.conda/pkgs" CONDA_ENVS_PATH="$CP"
export PATH="$HOME/envs/snakemake/bin:$PATH"
step(){ local n="$1" t="$2"; shift 2; [[ -f "$R/steps/$n.done" ]] && { echo "[skip] $n"; return 0; }
  echo "[step:$n] $(date +%H:%M:%S) start"; date +%s > "$R/steps/$n.start"
  timeout -k 60 "$t" "$@"; local rc=$?; [[ $rc -ne 0 ]] && { echo "[step:$n] FAILED rc=$rc"; return $rc; }
  date +%s > "$R/steps/$n.done"; echo "[step:$n] done"; }
cd "$PDIR" || exit 2
echo "=== snakemake 单株 snippy $(date) ==="
step snippy1 3600 bash -c "$SP --use-conda --conda-prefix $CP --cores 4 --rerun-incomplete results/snippy/GCF_000009065.1/snps.vcf"
echo "=== header 核验 ==="
step header 120 bash -c '
  A=results/snippy/GCF_000009065.1/snps.aligned.fa
  V=results/snippy/GCF_000009065.1/snps.vcf
  echo "--- aligned.fa headers ---"; grep "^>" "$A" | head -10
  echo "--- vcf CHROM ---"; cut -f1 "$V" | grep -v "^#" | sort -u
  echo "--- aligned.fa 长度 ---"; wc -c "$A"
'
echo "=== done $(date) ==="
