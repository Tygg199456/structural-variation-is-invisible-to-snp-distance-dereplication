# Zenodo metadata — the values to paste into the record

**Upload type**　Software

**Title**
Structural variation is invisible to SNP-distance dereplication: the reference strain
of *Yersinia pestis* merges with its own attenuated derivative

**Creators** (one row per author, in manuscript order; add affiliation and ORCID)
1. `<Family name>, <Given name>` — `<Affiliation>` — ORCID `<0000-0000-0000-0000>`
2. …

**Description** (plain text is fine; the 3–5 sentences below can be pasted as is)

Analysis code, derived tables and figure scripts for the manuscript *Structural
variation is invisible to SNP-distance dereplication: the reference strain of
Yersinia pestis merges with its own attenuated derivative*. The study maps public
*Y. pestis* assemblies to a single reference and shows that three measurements built
this way fail without raising an error: reference-based quality metrics act as
divergence filters, a SNP distance over jointly called sites cannot see a 102.6 kb
deletion, and absence relative to the reference has unknown polarity. The archive
contains the Snakemake pipeline, the six figure scripts, every derived table, a
per-number registry of the Results section (`results/numbers.json` plus
`docs/NUMBER_REGISTRY.md`), the supplementary tables, and the second-reference
replication. Raw assemblies and per-sample variant calls are not archived; they
regenerate from `pipeline/config.yaml` and the accessions listed in Table S2.

**Keywords**
Yersinia pestis; comparative genomics; dereplication; structural variation; reference
bias; reproducible analysis

**Licence**　Creative Commons Attribution 4.0 International (CC BY 4.0)
> the record contains derived data and figures; the code is MIT, and the split is
> stated in the repository README and in `LICENSE` / `LICENSE-DATA`.

**Version**　`<tag of the release>`

**Language**　English

**Related identifiers**
- `is identical to` → `https://github.com/Tygg199456/structural-variation-is-invisible-to-snp-distance-dereplication/tree/<tag>`
- `is documented by` → `https://github.com/Tygg199456/structural-variation-is-invisible-to-snp-distance-dereplication/blob/<tag>/README.md`

## Which DOI goes in the manuscript

Paste the **version DOI** of the published record, or the **concept DOI** (which
resolves to the latest version) — decide once and keep it consistent with the
repository README badge. `results/numbers.json` and the manuscript must agree.

The concept DOI also goes in the repository README, as the entry point to every
version.
