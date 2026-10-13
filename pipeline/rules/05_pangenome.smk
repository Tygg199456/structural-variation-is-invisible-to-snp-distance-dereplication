# ------------------------------------------------------------ 05_pangenome ---
# Tier A only. Draft assemblies are excluded by construction: a contig gap
# is indistinguishable from a gene deletion, and the whole study measures
# gene deletion. Fig.S3 quantifies exactly that bias.

rule panaroo:
    input:
        gffs = lambda w: expand(f"{RES}/annotate/{{acc}}/{{acc}}.gff3", acc=tier_a())
    output:
        rtab = f"{RES}/pangenome/gene_presence_absence.Rtab",
        csv  = f"{RES}/pangenome/gene_presence_absence.csv",
        aln  = f"{RES}/pangenome/core_gene_alignment.aln"
    threads: 32
    resources:
        mem_mb = 48000
    params:
        mode  = config["pangenome"]["clean_mode"],
        extra = config["pangenome"]["extra"],
        core  = config["pangenome"]["core_threshold"]
    log: f"{LOGS}/05_panaroo.log"
    conda: "../envs/panaroo.yaml"
    shell:
        r"""
        panaroo -i {input.gffs} -o {RES}/pangenome \
          --clean-mode {params.mode} {params.extra} \
          --core_threshold {params.core} -a core -t {threads} &> {log}
        """


rule ppanggolin:
    """Second opinion on persistent / shell / cloud partitioning."""
    input:
        gffs = lambda w: expand(f"{RES}/annotate/{{acc}}/{{acc}}.gff3", acc=tier_a())
    output:
        h5 = f"{RES}/pangenome/ppanggolin/pangenome.h5"
    threads: 32
    log: f"{LOGS}/05_ppanggolin.log"
    conda: "../envs/ppanggolin.yaml"
    shell:
        r"""
        ls {input.gffs} | awk '{{n=$0; sub(/.*\//,"",n); sub(/\.gff3$/,"",n);
          print n"\t"$0}}' > {RES}/pangenome/ppanggolin_input.tsv
        ppanggolin workflow --anno {RES}/pangenome/ppanggolin_input.tsv \
          -o {RES}/pangenome/ppanggolin -c {threads} --force &> {log}
        """


rule figS3_draft_pav_bias:
    """
    Fig.S3 -- the empirical justification for barring drafts from PAV.
    Take strains that have BOTH a complete and a draft assembly in NCBI,
    run each through the same panaroo settings, and count the losses the
    draft invents. This is the evidence, not an assertion.
    """
    input:
        meta = f"{RES}/qc/assembly_metadata.tsv",
        rtab = f"{RES}/pangenome/gene_presence_absence.Rtab"
    output:
        pdf = f"{RES}/figures/figS3_draft_pav_bias.pdf",
        tsv = f"{RES}/tables/figS3_spurious_losses.tsv"
    params:
        n_pairs = config["analysis"]["figS3_draft_pav_bias"]["n_paired_strains"]
    threads: 16
    log: f"{LOGS}/05_figS3.log"
    conda: "../envs/py.yaml"
    script: "../scripts/draft_pav_bias.py"
