# -------------------------------------------------------------- 10_figures ---
# Fig.1 carries a deliberate design decision: only populations that passed all
# four clock gates get a calendar axis. The rest are drawn with relative
# branch lengths on a greyed background. Methodological honesty becomes a
# visual element rather than a buried caveat.

rule fig1_tree:
    input:
        tree     = f"{RES}/dating/beast_global.mcc.tree",
        gates    = f"{RES}/dating/population_clock_gates.tsv",
        pav      = f"{RES}/pangenome/gene_presence_absence.Rtab",
        events   = f"{RES}/asr/pastml_events.tsv",
        clusters = f"{RES}/dedup/tierB_clusters.tsv",
        meta     = f"{RES}/qc/assembly_metadata.tsv"
    output:
        pdf = f"{RES}/figures/fig1_tree.pdf"
    log: f"{LOGS}/10_fig1.log"
    conda: "../envs/r.yaml"
    script: "../scripts/fig1_tree.R"


rule fig2_pangenome:
    input:
        rtab   = f"{RES}/pangenome/gene_presence_absence.Rtab",
        pseudo = f"{RES}/pseudogene/pseudo_counts.tsv",
        h5     = f"{RES}/pangenome/ppanggolin/pangenome.h5",
        meta   = f"{RES}/qc/assembly_metadata.tsv"
    output:
        pdf = f"{RES}/figures/fig2_pangenome.pdf"
    log: f"{LOGS}/10_fig2.log"
    conda: "../envs/py.yaml"
    script: "../scripts/fig2_pangenome.py"


rule fig3_parallel_loss:
    input:
        events = f"{RES}/asr/pastml_events.tsv",
        table2 = f"{RES}/tables/table2_parallel_genes.tsv"
    output:
        pdf = f"{RES}/figures/fig3_parallel_loss.pdf"
    log: f"{LOGS}/10_fig3.log"
    conda: "../envs/py.yaml"
    script: "../scripts/fig3_parallel_loss.py"


rule fig4_circos_is:
    input:
        overlap = f"{RES}/is_elements/breakpoint_overlap.tsv",
        ref     = f"{DATA}/ref/{REF}.fna"
    output:
        pdf = f"{RES}/figures/fig4_circos_is.pdf"
    log: f"{LOGS}/10_fig4.log"
    conda: "../envs/py.yaml"
    script: "../scripts/fig4_circos.py"


rule fig5_fba:
    input:
        growth = f"{RES}/fba/growth_matrix.tsv",
        auxo   = f"{RES}/fba/auxotrophy_matrix.tsv",
        gates  = f"{RES}/dating/population_clock_gates.tsv"
    output:
        pdf = f"{RES}/figures/fig5_fba.pdf"
    log: f"{LOGS}/10_fig5.log"
    conda: "../envs/py.yaml"
    script: "../scripts/fig5_fba.py"


rule fig6_decay_rate:
    """
    The direct answer to Q1: step change at speciation, or ongoing decay?
    Loss events per branch against node age. Nodes in populations that failed
    the clock gates are plotted on relative time and visually flagged.
    """
    input:
        events = f"{RES}/asr/pastml_events.tsv",
        tree   = f"{RES}/dating/beast_global.mcc.tree",
        gates  = f"{RES}/dating/population_clock_gates.tsv"
    output:
        pdf = f"{RES}/figures/fig6_decay_rate.pdf",
        tsv = f"{RES}/tables/fig6_regression.tsv"
    log: f"{LOGS}/10_fig6.log"
    conda: "../envs/py.yaml"
    script: "../scripts/fig6_decay_rate.py"


rule figS4_dedup_topology:
    input:
        full     = f"{RES}/phylogeny/allB_veryfasttree.nwk",
        backbone = f"{RES}/phylogeny/backbone.treefile"
    output:
        pdf = f"{RES}/figures/figS4_dedup_topology.pdf"
    log: f"{LOGS}/10_figS4.log"
    conda: "../envs/py.yaml"
    script: "../scripts/figS4_topology_compare.py"


rule table1_datasets:
    input:
        meta  = f"{RES}/qc/assembly_metadata.tsv",
        tierA = f"{RES}/dedup/tierA_keep.tsv",
        tierB = f"{RES}/dedup/tierB_backbone.tsv",
        tierC = f"{RES}/ancient/tierC_keep.tsv",
        qc    = f"{RES}/qc/tier_report.tsv"
    output:
        tsv = f"{RES}/tables/table1_datasets.tsv"
    log: f"{LOGS}/10_table1.log"
    conda: "../envs/py.yaml"
    script: "../scripts/table1_datasets.py"
