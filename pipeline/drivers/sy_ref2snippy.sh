#!/usr/bin/env bash
# sy_ref2snippy —— second-reference sensitivity run (v2, driven by snakemake so
# the per-rule conda env is activated properly).
#
# Maps the SAME Tier A set (results/qc/tierA_refined.tsv) against a second
# reference (GCF_000013825.1 = Y. pestis Antiqua, biovar Antiqua) and caches
# snippy output under results_ref2/.  Answers: does the reference choice, rather
# than the pipeline, drive the dereplication result?
set -uo pipefail
TASK="${1:-sy_ref2snippy}"
TASKDIR="$HOME/work/SY"; R="$TASKDIR/run/$TASK"
PDIR="$TASKDIR/pipeline/plague-decay"
REF=GCF_000013825.1
SP="$HOME/envs/snakemake/bin/snakemake"
CP="$HOME/envs/snakemake-envs"
DATASETS=<ENV_PREFIX>/snakemake-envs/cc3615509776b5f8fc7a3d2bb7e6e1ad_/bin/datasets
T="$TASKDIR/scratch/ref2"
mkdir -p "$TASKDIR/logs" "$R/steps" "$T"
date +%s > "$R/start"
trap 'echo "$? $(date +%s)" > "$R/exit"' EXIT
export PYTHONUNBUFFERED=1
export CONDA_PKGS_DIRS="$HOME/.conda/pkgs"
export CONDA_ENVS_PATH="$CP"

step(){ local n="$1" t="$2"; shift 2
  [[ -f "$R/steps/$n.done" ]] && { echo "[skip] $n"; return 0; }
  echo "[step:$n] $(date +%H:%M:%S) start"
  date +%s > "$R/steps/$n.start"
  timeout -k 60 "$t" "$@"; local rc=$?
  [[ $rc -ne 0 ]] && { echo "[step:$n] FAILED rc=$rc"; return $rc; }
  date +%s > "$R/steps/$n.done"; echo "[step:$n] $(date +%H:%M:%S) done"; }

cd "$PDIR" || { echo "CDFAIL $PDIR"; exit 2; }
echo "=== sy_ref2snippy v2 start $(date) ref=$REF ==="

step ref 1800 bash -c "set -euo pipefail
  if [ -s data/ref/$REF.fna ] && [ -s data/ref/$REF.gbff ]; then echo 'ref present'; exit 0; fi
  rm -rf $T/x $T/ref.zip
  $DATASETS download genome accession $REF --include genome,gbff --filename $T/ref.zip
  unzip -o -q $T/ref.zip -d $T/x
  cat $T/x/ncbi_dataset/data/$REF/*.fna  > data/ref/$REF.fna
  cat $T/x/ncbi_dataset/data/$REF/*.gbff > data/ref/$REF.gbff
  ls -la data/ref/$REF.fna data/ref/$REF.gbff"

step cfg 60 bash -c "set -euo pipefail
  sed -e '0,/^  accession: GCF_000009065.1/s//  accession: $REF/' \
      -e 's|^  results: results\$|  results: results_ref2|' \
      -e 's|^  logs: logs\$|  logs: logs_ref2|' \
      config.yaml > config.ref2.yaml
  grep -nE '^  (accession|results|logs):' config.ref2.yaml | head"

step targets 60 bash -c "set -euo pipefail
  cut -f1 results/qc/tierA_refined.tsv | tail -n +2 | grep -E '^GC[AF]_' | sort -u \
    | sed 's|^|results_ref2/snippy/|; s|\$|/snps.vcf|' > $T/targets.txt
  echo \"targets=\$(wc -l < $T/targets.txt)\""

step snippy 28800 bash -c "cd $PDIR
  $SP --configfile config.ref2.yaml --use-conda --conda-prefix $CP \
      --cores 24 --rerun-incomplete --latency-wait 60 --keep-going \
      \$(cat $T/targets.txt)"

step report 180 bash -c "
  echo '=== ref2 snippy outputs ==='
  ls -d $TASKDIR/results_ref2/snippy/*/ 2>/dev/null | wc -l
  ls -l $TASKDIR/results_ref2/snippy/GCF_000009065.1/snps.vcf $TASKDIR/results_ref2/snippy/GCF_001293415.1/snps.vcf 2>/dev/null || echo '(one of the CO92 pair missing)'
  echo '=== snakemake tail ==='
  tail -20 $TASKDIR/logs/sy_ref2snippy.out 2>/dev/null | tail -20"
echo "=== sy_ref2snippy v2 done $(date) ==="
