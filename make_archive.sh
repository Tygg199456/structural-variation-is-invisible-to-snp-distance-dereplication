#!/usr/bin/env bash
# Rebuild this archive from its two sources.  The repository IS the archive, so this
# script is documentation of how it was assembled (and how to refresh it), not a
# required step for a reader.
#
#   HOST   compute host holding the pipeline (ssh alias)
#   LOCAL  the working project directory holding the figure scripts and the result mirrors
#
#   HOST=<compute-host> LOCAL=~/<project-dir> ./make_archive.sh /tmp/archive
#
# It deliberately copies only code + derived tables.  It never touches
# results/snippy/, data/genomes/, db/ or any credential file.
set -euo pipefail
HOST="${HOST:-<compute-host>}"
LOCAL="${LOCAL:-$HOME/<project-dir>}"
OUT="${1:-$PWD/archive-rebuild}"

PIPELINE_FILES=(Snakefile config.yaml README.md rules scripts envs)
DRIVER_FILES=(sy_ref2snippy.sh sy_snippy1.sh sy_gate1c.sh sy_gate1d.sh sy_dedup_modes.sh sy_prof.sh)
RSYNC_EXCLUDES=(--exclude='__pycache__' --exclude='*.pyc' --exclude='*.bak*' --exclude='*.tab')

mkdir -p "$OUT"/{pipeline/drivers,figures/scripts,results/qc,results/dedup,results/figures,results/ref2,docs,env,tables,scripts}

echo "== pipeline code from $HOST =="
tar -C "$OUT/pipeline" -xf - <<<"$(
  ssh "$HOST" "cd ~/work/SY && tar chf - ${RSYNC_EXCLUDES[*]} \
      $(for f in "${PIPELINE_FILES[@]}"; do echo -n "pipeline/plague-decay/$f "; done)"
)"
# flatten pipeline/plague-decay/* -> pipeline/*
mv "$OUT/pipeline/plague-decay"/* "$OUT/pipeline/" && rmdir "$OUT/pipeline/plague-decay"
ssh "$HOST" "cd ~/work/SY/scripts && tar chf - ${DRIVER_FILES[*]}" | tar -C "$OUT/pipeline/drivers" -xf -

echo "== figure scripts and renders from $LOCAL =="
cp "$LOCAL"/scripts/fig/*.R "$LOCAL"/scripts/fig/prep_figure_data.py "$OUT/figures/scripts/"
cp "$LOCAL"/scripts/audit_numbers.py "$LOCAL"/scripts/_numbers.py "$OUT/scripts/"
cp "$LOCAL"/figures/SHA256SUMS.txt "$OUT/figures/"
for f in fig1_pipeline_schematic fig2_noise_floor fig3_quality_metrics \
         fig4_spatial_morphology fig5_structural_blindness fig6_polarity; do
  cp "$LOCAL/figures/$f.pdf" "$LOCAL/figures/$f.png" "$OUT/figures/"
done

echo "== derived tables =="
cp "$LOCAL"/data/qc/*.tsv "$LOCAL"/data/qc/*.txt "$OUT/results/qc/"
rm -f "$OUT/results/qc/METHODS_TODO.md" "$OUT/results/qc/_replinv_seqs.tsv"
cp "$LOCAL"/data/dedup/sites_tierA.tsv "$LOCAL"/data/dedup/sites_tierB.tsv \
   "$LOCAL"/data/dedup/snp_distance_raw_tierA.npy "$LOCAL"/data/dedup/snp_distance_raw_tierB.npy \
   "$LOCAL"/data/dedup/snp_distance_raw_tier*.order.txt "$LOCAL"/data/dedup/genotypes_tierA.npz \
   "$LOCAL"/data/dedup/tierA_clusters.*.tsv "$LOCAL"/data/dedup/tierB_backbone.tsv \
   "$LOCAL"/data/dedup/tierB_clusters.tsv "$OUT/results/dedup/"
cp "$LOCAL"/results/figures/*.csv "$LOCAL"/results/figures/*.prov "$OUT/results/figures/"
cp "$LOCAL"/results/numbers.json "$OUT/results/"
ssh "$HOST" 'cd ~/work/SY/pipeline/plague-decay/results_ref2 && tar chf - numbers.json verify.json' \
  | tar -C "$OUT/results/ref2" -xf -

echo "== scrub =="
bash "$(dirname "$0")/scripts/scrub_paths.sh" "$OUT"

echo "== manifest =="
( cd "$OUT" && find . -type f -not -name MANIFEST.sha256 -not -name .DS_Store | sort \
    | while read -r f; do printf '%s  %s\n' "$(shasum -a 256 "$f" | cut -d' ' -f1)" "${f#./}"; done > MANIFEST.sha256 )
echo "archive rebuilt at $OUT  ($(find "$OUT" -type f | wc -l) files)"
