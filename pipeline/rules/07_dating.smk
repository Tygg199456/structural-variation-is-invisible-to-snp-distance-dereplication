# --------------------------------------------------------------- 07_dating ---
# Temporal calibration.
#
# THE FAILURE MODE THAT DESTROYS THIS ANALYSIS SILENTLY:
# BEAST2 treats the supplied alignment as the complete set of sites. Feeding
# ~25,000 variable sites with no constant-site weights rescales branch lengths
# by roughly 4.65e6 / 2.5e4 ~= 190x. The run converges, Tracer looks healthy,
# and every rate and node age is wrong by two orders of magnitude.
# constantSiteWeights (AscertainedAlignment) is therefore mandatory and is
# asserted before the XML is written.
#
# Metadata support is adequate (measured): collection_date 89.1% on Tier A,
# 94.7% overall; year span 1939-2025 and 1900-2025. Both clear the 60% gate,
# so the full plan runs without degradation.

rule tempest_regression:
    """Root-to-tip regression per population; first of the four clock gates."""
    input:
        tree = f"{RES}/phylogeny/backbone.treefile",
        meta = f"{RES}/qc/assembly_metadata.tsv",
        pops = f"{RES}/phylogeny/population_assignment.tsv"
    output:
        tsv = f"{RES}/dating/tempest_by_population.tsv",
        pdf = f"{RES}/figures/figS5_root_to_tip.pdf"
    log: f"{LOGS}/07_tempest.log"
    conda: "../envs/py.yaml"
    script: "../scripts/root_to_tip.py"


rule population_clock_gates:
    """
    A population earns its own clock only if ALL four gates pass:
      n(dated tips) >= 30, year span >= 50, TempEst R2 >= 0.30 with positive
      slope, and a passing date-randomization test.
    Expected outcome: 1.ORI passes; 2.MED is borderline; the ancient and
    Pestoides lineages do not. Those get LSD2 relative dating and the
    manuscript reports no absolute ages for them.
    """
    input:
        tempest = f"{RES}/dating/tempest_by_population.tsv",
        drt     = f"{RES}/dating/drt_summary.tsv"
    output:
        tsv = f"{RES}/dating/population_clock_gates.tsv"
    params:
        gates = config["dating"]["population_clock_gates"]
    log: f"{LOGS}/07_clock_gates.log"
    conda: "../envs/py.yaml"
    script: "../scripts/clock_gates.py"


rule select_dating_taxa:
    """
    350 taxa (hard cap 400): ~60 ancient + ~170 dated Tier A + ~120 backbone
    supplements, the last chosen to maximise sampling-year span rather than
    assembly quality. Dating needs temporal spread; PAV needs Tier A. Keep
    those two selection criteria separate.
    """
    input:
        backbone = f"{RES}/dedup/tierB_backbone.tsv",
        tierA    = f"{RES}/dedup/tierA_keep.tsv",
        ancient  = f"{RES}/ancient/tierC_keep.tsv",
        meta     = f"{RES}/qc/assembly_metadata.tsv",
        hc5      = f"{RES}/phylogeny/hc5_clusters.tsv"
    output:
        tsv = f"{RES}/dating/dating_taxa.tsv"
    params:
        target      = config["dating"]["taxa_target"],
        cap         = config["dating"]["taxa_hard_cap"],
        composition = config["dating"]["composition"]
    log: f"{LOGS}/07_select_taxa.log"
    conda: "../envs/py.yaml"
    script: "../scripts/select_dating_taxa.py"


rule beast_xml:
    input:
        var    = f"{RES}/align/dating.var.aln",
        fconst = f"{RES}/align/dating.fconst.txt",
        taxa   = f"{RES}/dating/dating_taxa.tsv",
        dates  = f"{RES}/dating/tip_dates.tsv"
    output:
        xml = f"{RES}/dating/beast_global_{{chain}}.xml"
    params:
        clock  = config["dating"]["global"]["clock"],
        site   = config["dating"]["global"]["site_model"],
        prior  = config["dating"]["global"]["tree_prior"],
        length = config["dating"]["global"]["chain_length"],
        store  = config["dating"]["global"]["store_every"]
    log: f"{LOGS}/07_beast_xml_{{chain}}.log"
    conda: "../envs/py.yaml"
    script: "../scripts/make_beast_xml.py"   # asserts constantSiteWeights


rule beast_run:
    input:
        xml = f"{RES}/dating/beast_global_{{chain}}.xml"
    output:
        trees = f"{RES}/dating/beast_global_{{chain}}.trees",
        log   = f"{RES}/dating/beast_global_{{chain}}.log"
    threads: 8
    resources:
        runtime_h = 192          # 5-8 days CPU with BEAGLE SSE; 1-2 d on GPU
    params:
        beagle = config["dating"]["global"]["beagle"]
    log: f"{LOGS}/07_beast_run_{{chain}}.log"
    conda: "../envs/beast.yaml"
    shell:
        r"""
        beast -threads {threads} {params.beagle} -overwrite \
          -prefix {RES}/dating/ {input.xml} &> {log}
        """


rule beast_combine:
    input:
        trees = expand(f"{RES}/dating/beast_global_{{chain}}.trees", chain=[1, 2]),
        logs  = expand(f"{RES}/dating/beast_global_{{chain}}.log",   chain=[1, 2])
    output:
        mcc = f"{RES}/dating/beast_global.mcc.tree",
        ess = f"{RES}/dating/ess_summary.tsv"
    params:
        ess_min = config["dating"]["global"]["ess_threshold"]
    log: f"{LOGS}/07_beast_combine.log"
    conda: "../envs/beast.yaml"
    shell:
        r"""
        logcombiner -log {input.trees} -o {RES}/dating/combined.trees -burnin 10 &> {log}
        treeannotator -heights median -burnin 0 \
          {RES}/dating/combined.trees {output.mcc} &>> {log}
        python scripts/check_ess.py {input.logs} {params.ess_min} > {output.ess}
        """


rule date_randomization:
    """>=20 year-shuffled replicates; the 95% HPDs must not overlap the real run."""
    input:
        xml = f"{RES}/dating/beast_global_1.xml"
    output:
        tsv = f"{RES}/dating/drt_summary.tsv"
    params:
        reps = config["dating"]["population_clock_gates"]["date_randomization"]["replicates"]
    threads: 16
    log: f"{LOGS}/07_drt.log"
    conda: "../envs/beast.yaml"
    script: "../scripts/date_randomization.py"


rule lsd2_relative:
    """
    Fallback dating for populations that fail the gates, and the MVP path.
    Absolute ages are NOT reported for these; a fixed-rate scan documents the
    sensitivity of the relative ordering.
    """
    input:
        tree   = f"{RES}/phylogeny/{{set}}.treefile",
        dates  = f"{RES}/dating/tip_dates.tsv",
        fconst = f"{RES}/align/{{set}}.fconst.txt"
    output:
        nex = f"{RES}/dating/lsd2_{{set}}.nex"
    params:
        rates = config["dating"]["fallback"]["rate_scan"]
    log: f"{LOGS}/07_lsd2_{{set}}.log"
    conda: "../envs/iqtree.yaml"
    shell:
        r"""
        FC=$(cat {input.fconst})
        TOT=$(python -c "print(sum(int(x) for x in '$FC'.split(',')))")
        iqtree2 -s {RES}/align/{wildcards.set}.var.aln \
          -te {input.tree} --date {input.dates} --date-ci 100 \
          --dating LSD -fconst $FC \
          --prefix {RES}/dating/lsd2_{wildcards.set} -redo &> {log}
        """
