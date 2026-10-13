# Reproducing the figures and the numbers from this archive

Two paths, very different costs. Read the README's "Reproduction boundary" first.

## Path 1 — figures and numbers from the archived tables (minutes, no downloads)

```bash
cd <this repository>
python3 scripts/audit_numbers.py          # re-derives every Results number -> results/numbers.json
python3 scripts/make_number_registry.py   # renders docs/NUMBER_REGISTRY.md from it
python3 figures/scripts/prep_figure_data.py   # rebuilds results/figures/*.csv from results/qc, results/dedup
Rscript figures/scripts/fig1_pipeline_schematic.R   # ... and fig2 … fig6 likewise
```

Environment for path 1:
- Python 3 with `numpy` (the audit and the prep script also need nothing else).
- R with `ggplot2`, `patchwork`, `ggrepel`, `scales`, `ragg`, `data.table`
  (`env/R_SESSION_INFO.txt` records the exact versions used for the published figures).

`figures/scripts/prep_figure_data.py` writes each CSV by `os.replace`, so a failed step
cannot leave a half file for R to read, and raises rather than skipping.

## Path 2 — derived tables from public accessions (the real cost)

```bash
snakemake --use-conda -c 32 --configfile config.yaml gate1c   # Tier A: download, QC gate, snippy, matrix, dedup
```

Needs NCBI Datasets CLI, ~157 GB of scratch, a bakta database (39,401,553,410 bytes,
MD5 `f8823533b789dd315025fdcc46f1c1…` as recorded in `config.yaml`) and several CPU-days.
The per-sample `snippy` output and the raw assemblies are not archived; they regenerate
from `config.yaml` and the accessions in the manuscript's Table S2.

The second-reference replication (`results/ref2/`) was produced by
`pipeline/drivers/sy_ref2snippy.sh` followed by
`python3 pipeline/scripts/ref2_analysis.py`; it maps the same 256 genomes against
*Y. pestis* Antiqua (`GCF_000013825.1`) and re-applies the same threshold rule on the
same sweep grid.

## Assertions that are meant to fail loudly

- `pipeline/scripts/genotype_matrix.py` exits non-zero if any sample's
  `snps.aligned.fa` is missing — a VCF alone cannot distinguish REF from unknown.
- `pipeline/scripts/snp_distance.py` applies the joint-call floor before any threshold.
- `figures/scripts/fig2_noise_floor.R` stops if the CO92-against-itself alignment
  produces an uncalled run of 10 kb or more, which would invalidate the spatial test.
- `figures/scripts/fig1_pipeline_schematic.R` stops if any illustrative site is
  `called_diff` in both genomes, a state the three-state cartoon cannot encode.
- `scripts/audit_numbers.py` raises on a missing input rather than skipping it.
