#!/usr/bin/env python3
"""
Create contract stubs for the scripts not yet written out in full.

Each stub records its input/output contract and the one analytical decision
that matters for it, so the implementation cannot drift from the design.
Stubs raise NotImplementedError rather than returning silently wrong results.
Run once: python scripts/_make_stubs.py
"""
import pathlib

HERE = pathlib.Path(__file__).resolve().parent

# name -> (one-line purpose, the decision that must not be lost)
STUBS = {
    "flatten_metadata.py": (
        "NCBI Datasets JSONL -> flat TSV + field completeness report.",
        "Normalise missing tokens from config (metadata.missing_tokens). "
        "Measured targets to reproduce: Tier A collection_date 89.1%, "
        "geo_loc_name 94.2%, host 60.1%, isolation_source 39.1%; "
        "all-set 94.7 / 94.0 / 29.4 / 66.6%."),
    "genotype_matrix.py": (
        "Per-sample snippy VCFs -> int8 genotype matrix over the union of sites.",
        "Code 0=missing, 1..4=A,C,G,T. Never build a whole-genome alignment "
        "here: 3453 x 4.65 Mb is ~16 GB, this matrix is ~86 MB."),
    "figS1_checkm_bias.py": (
        "Fig.S1 -- CheckM completeness vs pseudogene count, with regression.",
        "A significant negative slope is the published justification for not "
        "filtering on completeness. Report slope, R2 and p in the table."),
    "tier_c_select.py": (
        "Select ancient libraries passing depth + damage authentication.",
        "UDG-treated libraries carry little terminal deamination and will "
        "fail the damage gate; read library prep from AncientMetagenomeDir "
        "and exempt them rather than silently dropping real data."),
    "coarse_pav.py": (
        "Breadth-based gene presence for ancient genomes.",
        "Presence = breadth >= 90% of gene length at >= 3x. Emit presence/"
        "absence/UNKNOWN; never emit a frameshift or pseudogene call."),
    "downsample_control.py": (
        "Fig.S2 -- false-loss rate when Tier A is downsampled to Tier C depth.",
        "Without this control the ancient losses in Fig.3 are indefensible. "
        "Match the empirical Tier C depth distribution, not a single value."),
    "pseudo_matrix.py": (
        "Merge pseudofinder + bakta pseudo calls into a strain x gene matrix.",
        "Report only calls supported by BOTH tools. Lower recall, far lower "
        "false positive rate -- correct trade when the count IS the result."),
    "is_breakpoints.py": (
        "Co-localisation of IS copies with loss/pseudogenisation breakpoints.",
        "Test against a position-randomised null; raw overlap counts mean "
        "nothing when IS elements occupy a large share of the genome."),
    "draft_pav_bias.py": (
        "Fig.S3 -- spurious losses introduced by draft assemblies.",
        "Use strains with BOTH complete and draft assemblies in NCBI; run "
        "each through identical panaroo settings. This is the evidence for "
        "barring Tier B from PAV, not an assertion."),
    "mask_alignment.py": (
        "Apply the combined mask BED to a full-length alignment.",
        "Write N, do not delete columns: snp-sites -C must still see the "
        "true constant-site counts for -fconst and constantSiteWeights."),
    "root_to_tip.py": (
        "Per-population root-to-tip regression (clock gate 3).",
        "Run per population, never on the whole species -- the global "
        "Y. pestis clock is unstable and a pooled regression hides it."),
    "clock_gates.py": (
        "Apply the four clock gates and decide BEAST vs LSD2 per population.",
        "All four must pass: n>=30 dated tips, span>=50 y, R2>=0.30 with "
        "positive slope, and a passing date-randomization test."),
    "select_dating_taxa.py": (
        "Assemble the ~350-taxon dating set.",
        "Backbone supplements are chosen to MAXIMISE sampling-year span, not "
        "assembly quality. Dating needs temporal spread; PAV needs Tier A."),
    "date_randomization.py": (
        "Date-randomization test, >=20 year-shuffled replicates.",
        "Pass = the real 95% HPD does not overlap the shuffled distribution."),
    "check_ess.py": (
        "Parse BEAST logs, assert ESS >= threshold for every parameter.",
        "Exit non-zero on failure; a converged-looking run with low ESS on "
        "one parameter is still unusable."),
    "prune_tree.py": (
        "Prune the dated MCC tree to Tier A + Tier C for ASR.",
        "Keep dated branch lengths: loss events must land on a real time "
        "axis even though only Tier A carries trustworthy PAV."),
    "combine_pav.py": (
        "Merge Tier A assembly PAV, Tier C coarse PAV and pseudogene states.",
        "Three-state output (present / absent / unknown). Tier C unknowns "
        "must stay unknown, not be imputed to absent."),
    "extract_events.py": (
        "PastML work dir -> tidy gene/branch/event table.",
        "Count a loss once per branch; repeated marginal states on one branch "
        "are not independent events."),
    "fetch_ipc815.py": (
        "Download iPC815 SBML (Charusanti 2011) from BiGG or BioModels.",
        "Verify 815 genes / 678 proteins / 963 metabolites / 1678 reactions "
        "after load; fail loudly if the model does not match."),
    "define_media.py": (
        "Build the three in silico media as exchange-bound sets.",
        "Mammalian blood must set the iron exchange lower bound to 0 -- host "
        "iron sequestration is what makes ybt/HPI necessary."),
    "strain_fba.py": (
        "Per-strain and per-ancestral-node knockout + pFBA on iPC815.",
        "Let GPR rules propagate; do not knock out reactions directly. "
        "Record which lost genes have NO GPR -- that is the 20% coverage "
        "caveat made quantitative."),
    "model_crossval.py": (
        "gapseq / CarveMe reconstructions for 10 representatives.",
        "If de novo models contradict the iPC815 trend, demote the FBA panel "
        "to supplementary rather than arguing around it."),
    "fig2_pangenome.py": (
        "Fig.2 -- pangenome curve, partitions, pseudogene distribution.",
        "State in the caption that PAV and pseudogenes come from Tier A only."),
    "fig3_parallel_loss.py": (
        "Fig.3 -- PastML event map + parallel-loss Manhattan plot.", ""),
    "fig4_circos.py": (
        "Fig.4 -- circos of four replicons, IS sites and loss breakpoints.", ""),
    "fig5_fba.py": (
        "Fig.5 -- lineage x medium growth heatmap, auxotrophy onset.",
        "Mark lineages lacking a calibrated clock; their onset times are "
        "relative, not calendar dates."),
    "fig6_decay_rate.py": (
        "Fig.6 -- loss events per branch vs node age (step vs continuous).",
        "Plot gate-failing populations on relative time and flag them "
        "visually; do not mix calendar and relative axes silently."),
    "figS4_topology_compare.py": (
        "Fig.S4 -- all-3453 VeryFastTree vs deduplicated backbone.",
        "Report Robinson-Foulds distance and the support-value distributions "
        "side by side. This answers 'why discard 1800 genomes'."),
    "table1_datasets.py": (
        "Table 1 -- full dataset inventory with tier, QC and metadata.", ""),
}

R_STUBS = {
    "fig1_tree.R": (
        "Fig.1 -- dated tree with PAV ring, event colouring, cluster bubbles.",
        "Only populations passing all four clock gates get a calendar axis; "
        "the rest are drawn on relative branch lengths and greyed. Make the "
        "methodological limit a visual element, not a buried caveat."),
    "enrichment.R": (
        "KEGG (ypo) / GO / COG enrichment of significant parallel-loss genes.",
        "Background = the Tier A pangenome, not the whole KEGG universe."),
}

PY_TMPL = '''#!/usr/bin/env python3
"""
{purpose}

{decision_block}
Inputs / outputs are supplied by Snakemake via the `snakemake` object.
"""

snakemake = globals().get("snakemake")  # noqa: F821


def main():
    raise NotImplementedError(
        "{name}: not implemented. See the module docstring for the contract."
    )


if __name__ == "__main__":
    main()
'''

R_TMPL = '''#!/usr/bin/env Rscript
# {purpose}
#
{decision_block}
# Inputs / outputs come from snakemake@input / snakemake@output.

stop("{name}: not implemented. See the header for the contract.")
'''


def block(text, prefix=""):
    if not text:
        return ""
    import textwrap
    wrapped = textwrap.fill(text, 74)
    if prefix:
        wrapped = "\n".join(prefix + ln for ln in wrapped.splitlines())
    return "DECISION THAT MUST NOT BE LOST:\n" + wrapped + "\n" if not prefix \
        else "# DECISION THAT MUST NOT BE LOST:\n" + wrapped + "\n"


written = 0
for name, (purpose, decision) in STUBS.items():
    p = HERE / name
    if p.exists():
        continue
    p.write_text(PY_TMPL.format(purpose=purpose, name=name,
                                decision_block=block(decision)))
    written += 1

for name, (purpose, decision) in R_STUBS.items():
    p = HERE / name
    if p.exists():
        continue
    p.write_text(R_TMPL.format(purpose=purpose, name=name,
                               decision_block=block(decision, prefix="# ")))
    written += 1

print(f"wrote {written} stub scripts to {HERE}")
