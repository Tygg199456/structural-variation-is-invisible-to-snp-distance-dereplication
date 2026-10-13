# plague-decay

Reductive evolution in *Yersinia pestis*: phylogenomics, ancestral gene loss,
and genome-scale metabolic consequence.

All tier sizes below are derived from measured NCBI Datasets output
(datasets CLI 18.38.0, taxid 632, 2026-10-09), not estimates.

---

## The three decisions this pipeline encodes

### 1. CheckM completeness is not a quality filter here

110 of 257 Tier A candidates (42.8%) fall below 98% completeness, and the
strains that drop out are the deep-branching anchors the study depends on
(Angola 96.76, Antiqua 97.70, Pestoides F ~97–98).

CheckM scores a genome against Enterobacteriaceae single-copy markers.
*Y. pestis* has pseudogenised or deleted a sizeable share of exactly those
markers — that is the phenomenon under study. Filtering on completeness
deletes the signal and calls it QC.

- **Contamination is kept** (≤2%): it detects mixed/chimeric assemblies, a
  genuinely separate problem.
- **Completeness is kept as a covariate** and tested against pseudogene count
  in Fig.S1.
- **Assembly completeness is measured marker-free**: `dnadiff` AlignedBases
  ≥95% against the CO92 chromosome.

Enforced in `scripts/assign_tiers.py`.

### 2. Tiers have single permitted uses, enforced in code

| Tier | n (measured / projected) | Permitted | Forbidden |
|---|---|---|---|
| **A** decay set | 257 candidates → ~210 after SNP dedup | PAV, pseudogenes, IS, plasmids, FBA, Gubbins | — |
| **B** backbone | 3,453 after GCA/GCF dedup → 1,300–1,900 backbone | SNP phylogeny, dating tips | **PAV, pseudogenes, FBA** |
| **C** ancient | 30 pilot (≥10×) → 60–80 full (≥5×) | SNP, deep calibration, coarse PAV | **frameshift, pseudogenes** |

A contig gap is indistinguishable from a gene deletion, and gene deletion is
the measurement. Fig.S3 quantifies that bias empirically rather than asserting
it. `assert_use()` in the Snakefile raises on violation.

### 3. Constant sites are counted, not conditioned away

`snp-sites -C` gives the exact A,C,G,T constant-site counts, which feed both:

- **IQ-TREE `-fconst`**, not `+ASC`. `+ASC` conditions the likelihood on
  constant sites being *unobservable*; here they are observed and counted, so
  `+ASC` is the wrong model and inflates branch lengths.
- **BEAST2 `constantSiteWeights`**. This one is a silent killer: feeding
  ~25k variable sites with no weights inflates branch lengths ~186×. The run
  converges, Tracer looks healthy, and every rate and node age is wrong by two
  orders of magnitude. `scripts/make_beast_xml.py` asserts the weights are
  present and plausible, and refuses to write the XML otherwise.

---

## Running it

```bash
mamba install -c conda-forge -c bioconda snakemake'>=8'
snakemake -n                                   # dry run

snakemake --use-conda -c32 tier_counts         # GATE 1 — stop and read
snakemake --use-conda -c32 mvp                 # ~5–7 days
snakemake --use-conda -c32 all                 # ~3 weeks
```

**Do not skip Gate 1.** It writes the realised n per tier so the placeholders
in `config.yaml` (Tier A 210, backbone 1300–1900) can be confirmed or
corrected before any long-running step starts.

### Before the first run

1. `config.yaml` → `annotation.db_version`: pin the exact Bakta DB release.
   One DB version for every genome, references included. Mixing PGAP vintages
   is the largest source of spurious pseudogene differences, and it is
   invisible downstream.
2. `config.yaml` → `whitelist.*.population`: fill from the supplementary
   tables of Morelli 2010 and Cui 2013. Do not type these from memory.
3. `config.yaml` → `outgroup.extra`: verify the two accessions resolve.
4. BEAST2 packages: `packagemanager -add ORC -add RLC -add bModelTest`
5. `$BAKTA_DB` must be set.
6. Build the pseudofinder DIAMOND reference db (*Y. pseudotuberculosis* plus
   related Enterobacterales).

---

## Timings (32 cores / 128 GB)

| Step | Scale | Time | Disk |
|---|---|---|---|
| Download assemblies | 4,585 | 2–4 h | ~25 GB |
| Tier B QC | ~4,500 | 8–12 h | — |
| Bakta (Tier A) | ~210 | 1.5–2 h | 5 GB |
| Panaroo | ~210 | 2–4 h, 32 GB peak | 10 GB |
| Pseudofinder | ~210 | 4–6 h | 5 GB |
| snippy (Tier B, `--ctgs`) | ~3,453 | 4–6 h | 60 GB |
| Genotype matrix + distance | 3,453 × 25k | ~10 min | 86 MB |
| Gubbins (**Tier A only**) | 257 | 2–4 h | — |
| IQ-TREE backbone | ~1,600 taxa | 12–36 h | — |
| VeryFastTree (Fig.S4) | 3,453 | 1–2 h | — |
| nf-core/eager | 30 → 80 libs | 1.5–3 d | **0.2–2 TB** |
| BEAST ×2 chains | 350 taxa | 5–8 d CPU / 1–2 d GPU | — |
| PastML + null (1000 sims) | — | 4–8 h | — |
| FBA (iPC815) | 210 × 3 media | **<30 min** | — |

Ancient DNA reads are the only storage problem. Run the ≥10× pilot (~30
libraries, 200–400 GB) end to end before committing disk for phase 2.

---

## Layout

```
Snakefile              targets: tier_counts | mvp | all
config.yaml            every locked parameter, with the reasoning inline
rules/00..10_*.smk     fetch, QC, dedup, ancient, annotate, pangenome,
                       phylogeny, dating, ASR, FBA, figures
scripts/               analysis scripts (4 written in full, 31 contract stubs)
envs/                  22 conda environments
```

Fully implemented: `assign_tiers.py`, `snp_dedup.py`, `snp_distance.py`,
`make_beast_xml.py`, `parallel_loss_null.py`.
The remaining 31 are contract stubs — each carries its I/O contract and the
one analytical decision that must survive implementation, and raises
`NotImplementedError` rather than returning quietly wrong results.

---

## Caveats to carry into the manuscript

- **iPC815 covers ~815 of ~4,000 genes (~20%).** Many lost genes have no GPR
  and are invisible to FBA. Scope conclusions to metabolic *modules*, state
  the coverage in Methods, and cross-check with gapseq/CarveMe. If the de novo
  models contradict the iPC815 trend, demote the FBA panel to supplementary.
- **The global molecular clock is unstable.** Only populations clearing all
  four gates get calendar ages; the rest get LSD2 relative dating and no
  absolute ages are reported for them. Fig.1 encodes this visually.
- **DURC.** *Y. pestis* is a Tier 1 Select Agent. Avoid framing results as a
  route to enhanced virulence or vaccine escape; include a DURC statement;
  budget 2–6 extra weeks for biosecurity review at mBio / PLoS Pathogens /
  Nature-family journals.
