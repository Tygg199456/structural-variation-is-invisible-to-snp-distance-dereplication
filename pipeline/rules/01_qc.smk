# --------------------------------------------------- 01_qc.smk  (PATCH v2) ---
# REPLACES the v1 file entirely.
#
# v1 made Gate 1 the heaviest step in the pipeline: `assign_tiers` required a
# dnadiff run for every one of the 4,585 records, so "check the counts before
# committing" meant 12-20 h and 25 GB of downloads first. That was a design
# error, not an intended trade-off.
#
# v2 splits the gate:
#
#   Gate 1a  assign_tiers_metadata  -- metadata only, no downloads, minutes.
#            Emits Tier A (~257) and Tier B (~3,453) candidate counts. This
#            is the gate you actually read before committing compute.
#
#   Gate 1b  dnadiff_vs_ref         -- Tier A candidates ONLY (~257, 1-2 h).
#
#   Tier B coverage comes from snippy's own snps.aligned.fa, which is
#   produced anyway. ~3,200 dnadiff runs eliminated as pure redundancy.
#
# GCA/GCF pairing is resolved from the `paired_accession` metadata field
# BEFORE anything is downloaded, saving ~1,132 downloads.
#
# Unchanged from v1 and still the central decision: CheckM completeness is
# recorded as a covariate, never used as a filter. 110/257 Tier A candidates
# (42.8%) fall below 98%, and the strains that drop out are the deep-branching
# anchors the study depends on.

checkpoint assign_tiers_metadata:
    """
    GATE 1a -- metadata only. No downloads, no alignments.
    Criteria available without sequence: assembly level, total length,
    CheckM contamination, contig count, N50, ambiguous fraction, and
    GCA/GCF pairing.
    """
    input:
        meta = f"{RES}/qc/assembly_metadata.tsv"
    output:
        tierA  = f"{RES}/qc/tierA_candidates.tsv",
        tierB  = f"{RES}/qc/tierB_candidates.tsv",
        report = f"{RES}/qc/tier_report.tsv",
        excl   = f"{RES}/qc/excluded_with_reason.tsv"
    params:
        cfg = config
    log: f"{LOGS}/01_assign_tiers_metadata.log"
    conda: "../envs/py.yaml"
    script: "../scripts/assign_tiers.py"      # run with stage="metadata"


rule dnadiff_vs_ref:
    """
    GATE 1b -- marker-free completeness for TIER A ONLY.
    Fraction of the CO92 chromosome covered. Tier B does not come through
    here; its coverage is read off snps.aligned.fa instead.
    """
    input:
        qry = f"{DATA}/genomes/{{acc}}/{{acc}}.fna",
        ref = f"{DATA}/ref/{REF}.fna"
    output:
        rpt = f"{RES}/qc/dnadiff/{{acc}}.report"
    log: f"{LOGS}/01_dnadiff/{{acc}}.log"
    conda: "../envs/mummer.yaml"
    shell:
        r"""
        d=$(dirname {output.rpt}); mkdir -p $d
        dnadiff -p $d/{wildcards.acc} {input.ref} {input.qry} &> {log}
        """


checkpoint refine_tier_a:
    """Apply the dnadiff gate to the Tier A candidate list."""
    input:
        cand    = f"{RES}/qc/tierA_candidates.tsv",
        dnadiff = lambda w: expand(f"{RES}/qc/dnadiff/{{acc}}.report",
                                   acc=tier_a_candidates())
    output:
        keep = f"{RES}/qc/tierA_refined.tsv"
    params:
        cfg = config
    log: f"{LOGS}/01_refine_tier_a.log"
    conda: "../envs/py.yaml"
    script: "../scripts/refine_tier_a.py"


rule tier_b_coverage_filter:
    """
    Tier B coverage from snippy, not dnadiff.
    snps.aligned.fa is reference-length, so the called fraction is already
    computed by genotype_matrix.py and stored in genotypes.npz.
    """
    input:
        npz  = f"{RES}/dedup/genotypes_tierB.npz",
        cand = f"{RES}/qc/tierB_candidates.tsv"
    output:
        keep = f"{RES}/qc/tierB_coverage_pass.tsv"
    params:
        min_aligned = lambda w: config["tiers"]["B"]["qc"]["min_dnadiff_aligned_pct"] / 100.0,
        whitelist   = config["whitelist"]
    log: f"{LOGS}/01_tier_b_coverage.log"
    conda: "../envs/py.yaml"
    script: "../scripts/tier_b_coverage.py"


rule figS1_checkm_bias:
    """
    Fig.S1 -- CheckM completeness vs pseudogene count.
    A significant negative slope is the published justification for not
    filtering on completeness.
    """
    input:
        qc     = f"{RES}/qc/tier_report.tsv",
        pseudo = f"{RES}/pseudogene/pseudo_counts.tsv"
    output:
        pdf = f"{RES}/figures/figS1_checkm_bias.pdf",
        tsv = f"{RES}/tables/figS1_regression.tsv"
    conda: "../envs/py.yaml"
    script: "../scripts/figS1_checkm_bias.py"


# ---------------------------------------------------------------- helpers ---
def _col(path, name="accession"):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    if not rows:
        return []
    try:
        i = rows[0].index(name)
    except ValueError:
        i = 0
    return [r[i] for r in rows[1:]]

def tier_a_candidates():
    return _col(checkpoints.assign_tiers_metadata.get().output.tierA)

def tier_b_candidates():
    return _col(checkpoints.assign_tiers_metadata.get().output.tierB)
