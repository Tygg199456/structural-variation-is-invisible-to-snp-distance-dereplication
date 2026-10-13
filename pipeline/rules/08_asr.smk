# ------------------------------------------------------------------ 08_asr ---
# Ancestral reconstruction of gene loss, then the null model that separates
# adaptation from drift.
#
# The tree used here is the DATED tree pruned to Tier A + Tier C. Branch
# lengths and the time axis are inherited from the full dated analysis, so
# loss events land on a real time axis even though only Tier A carries
# trustworthy PAV.

rule prune_tree_for_asr:
    input:
        tree  = f"{RES}/dating/beast_global.mcc.tree",
        tierA = f"{RES}/dedup/tierA_keep.tsv",
        tierC = f"{RES}/ancient/tierC_keep.tsv"
    output:
        nwk = f"{RES}/asr/dated_pruned_AC.nwk"
    log: f"{LOGS}/08_prune.log"
    conda: "../envs/py.yaml"
    script: "../scripts/prune_tree.py"


rule combine_pav:
    """Tier A assembly-based PAV + Tier C breadth-based coarse PAV."""
    input:
        rtab   = f"{RES}/pangenome/gene_presence_absence.Rtab",
        csv    = f"{RES}/pangenome/gene_presence_absence.csv",
        coarse = f"{RES}/ancient/coarse_pav.tsv",
        pseudo = f"{RES}/pseudogene/pseudo_calls_long.tsv"
    output:
        tsv = f"{RES}/asr/combined_states.tsv"
    log: f"{LOGS}/08_combine_pav.log"
    conda: "../envs/py.yaml"
    script: "../scripts/combine_pav.py"


rule pastml:
    input:
        tree   = f"{RES}/asr/dated_pruned_AC.nwk",
        states = f"{RES}/asr/combined_states.tsv"
    output:
        events = f"{RES}/asr/pastml_events.tsv",
        html   = f"{RES}/asr/pastml_map.html"
    threads: 16
    params:
        method = config["asr"]["prediction_method"],
        model  = config["asr"]["model"]
    log: f"{LOGS}/08_pastml.log"
    conda: "../envs/pastml.yaml"
    shell:
        r"""
        pastml --tree {input.tree} --data {input.states} \
          --prediction_method {params.method} --model {params.model} \
          --html {output.html} --work_dir {RES}/asr/pastml_work \
          --threads {threads} &> {log}
        python scripts/extract_events.py {RES}/asr/pastml_work > {output.events}
        """


rule parallel_loss_null:
    """
    Per gene: count independent losses on the tree, then simulate 1000 times
    at the ML loss rate on the SAME tree. Genes lost more often than the null
    allows are candidate adaptive losses; the rest are drift.
    This test is what makes the parallelism claim falsifiable -- the existing
    literature asserts convergence without a null.
    """
    input:
        events = f"{RES}/asr/pastml_events.tsv",
        tree   = f"{RES}/asr/dated_pruned_AC.nwk"
    output:
        tsv = f"{RES}/tables/table2_parallel_genes.tsv"
    threads: 16
    params:
        sims  = config["parallel_loss"]["simulations"],
        alpha = config["parallel_loss"]["fdr_alpha"]
    log: f"{LOGS}/08_parallel_null.log"
    conda: "../envs/py.yaml"
    script: "../scripts/parallel_loss_null.py"


rule enrichment:
    input:
        genes = f"{RES}/tables/table2_parallel_genes.tsv",
        anno  = f"{RES}/annotate/{REF}/{REF}.tsv"
    output:
        tsv = f"{RES}/asr/enrichment.tsv"
    log: f"{LOGS}/08_enrichment.log"
    conda: "../envs/r.yaml"
    script: "../scripts/enrichment.R"
