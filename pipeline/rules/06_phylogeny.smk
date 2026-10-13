# ------------------------------------------------------------ 06_phylogeny ---
# Masking, then ML trees.
#
# Two corrections over the naive recipe:
#
# 1. Gubbins runs on Tier A ONLY. Gubbins needs a FULL-LENGTH alignment, not a
#    SNP alignment. Tier A is 257 x 4.65 Mb (~1.2 GB) and is fine. Tier B
#    would be 3453 x 4.65 Mb (~16 GB) and the internal RAxML iterations will
#    exhaust memory. The recombinant regions Gubbins finds on Tier A are
#    merged into the mask BED and propagated to Tier B.
#
# 2. -fconst, NOT +ASC. snippy-core preserves full reference coordinates, so
#    the constant-site counts are known exactly and are supplied directly.
#    +ASC conditions the likelihood on constant sites being unobservable,
#    which is the wrong model when they have in fact been observed and
#    counted; it biases branch lengths upward. +ASC also aborts outright if
#    any invariant column survives into the alignment.

rule snippy_core:
    input:
        vcfs = lambda w: expand(f"{RES}/snippy/{{acc}}/snps.vcf",
                                acc=(tier_a() if w.set == "tierA" else tier_b())),
        ref  = f"{DATA}/ref/{REF}.gbff"
    output:
        full = f"{RES}/align/{{set}}.full.aln",
        aln  = f"{RES}/align/{{set}}.aln"
    log: f"{LOGS}/06_snippy_core_{{set}}.log"
    conda: "../envs/snippy.yaml"
    shell:
        r"""
        dirs=$(dirname {input.vcfs} | tr '\n' ' ')
        snippy-core --ref {input.ref} \
          --prefix {RES}/align/{wildcards.set} $dirs &> {log}
        """


rule build_mask:
    """Prophage + self-self repeats + IS elements + Gubbins recombinant regions."""
    input:
        ref      = f"{DATA}/ref/{REF}.fna",
        phastest = f"{RES}/mask/phastest_regions.bed",
        repeats  = f"{RES}/mask/repeats.bed",
        is_bed   = f"{RES}/mask/is_elements.bed",
        gubbins  = f"{RES}/recomb/tierA.recombination_predictions.gff"
    output:
        bed = f"{RES}/mask/combined_mask.bed"
    log: f"{LOGS}/06_build_mask.log"
    conda: "../envs/bedtools.yaml"
    shell:
        r"""
        {{
          cat {input.phastest} {input.repeats} {input.is_bed}
          awk 'BEGIN{{OFS="\t"}} !/^#/ {{print $1,$4-1,$5}}' {input.gubbins}
        }} | sort -k1,1 -k2,2n | bedtools merge -i - > {output.bed} 2> {log}
        wc -l {output.bed} >> {log}
        """


rule gubbins_tier_a:
    input:
        aln = f"{RES}/align/tierA.full.aln"
    output:
        gff = f"{RES}/recomb/tierA.recombination_predictions.gff"
    threads: 32
    log: f"{LOGS}/06_gubbins_tierA.log"
    conda: "../envs/gubbins.yaml"
    shell:
        r"""
        # Tier A only -- see header note. Expected result in Y. pestis is
        # near-zero recombination; report it as a negative control either way.
        run_gubbins.py --threads {threads} \
          --prefix {RES}/recomb/tierA {input.aln} &> {log}
        """


rule mask_alignment:
    input:
        aln = f"{RES}/align/{{set}}.full.aln",
        bed = f"{RES}/mask/combined_mask.bed"
    output:
        aln = f"{RES}/align/{{set}}.masked.full.aln"
    log: f"{LOGS}/06_mask_{{set}}.log"
    conda: "../envs/py.yaml"
    script: "../scripts/mask_alignment.py"


rule constant_sites:
    """
    snp-sites -C gives the exact A,C,G,T constant-site counts.
    These feed BOTH iqtree -fconst and the BEAST2 constantSiteWeights.
    """
    input:
        aln = f"{RES}/align/{{set}}.masked.full.aln"
    output:
        fconst = f"{RES}/align/{{set}}.fconst.txt",
        var    = f"{RES}/align/{{set}}.var.aln"
    log: f"{LOGS}/06_constant_sites_{{set}}.log"
    conda: "../envs/snpsites.yaml"
    shell:
        r"""
        snp-sites -C {input.aln} > {output.fconst} 2> {log}
        snp-sites -c {input.aln} > {output.var}   2>> {log}
        echo "constant sites A,C,G,T = $(cat {output.fconst})" >> {log}
        echo "variable sites = $(grep -c '^>' {output.var} >/dev/null; \
              awk 'NR==2{{print length($0)}}' {output.var})" >> {log}
        """


rule iqtree:
    input:
        var    = f"{RES}/align/{{set}}.var.aln",
        fconst = f"{RES}/align/{{set}}.fconst.txt"
    output:
        tree = f"{RES}/phylogeny/{{set}}.treefile"
    threads: 32
    params:
        model = config["phylogeny"]["model"],
        ufb   = config["phylogeny"]["ufboot"],
        alrt  = config["phylogeny"]["alrt"]
    log: f"{LOGS}/06_iqtree_{{set}}.log"
    conda: "../envs/iqtree.yaml"
    shell:
        r"""
        FC=$(cat {input.fconst})
        iqtree2 -s {input.var} -m {params.model} -fconst $FC \
          -B {params.ufb} --alrt {params.alrt} -T {threads} \
          --prefix {RES}/phylogeny/{wildcards.set} -redo &> {log}
        """


rule veryfasttree_all_b:
    """
    Fig.S4 -- all 3453 Tier B taxa, undeduplicated, ~1-2 h.
    Shows the backbone topology is unchanged by dedup while support values on
    the full set collapse. This is the answer to 'why did you discard 1800
    genomes', and it costs almost nothing.
    """
    input:
        var = f"{RES}/align/tierB_full.var.aln"
    output:
        nwk = f"{RES}/phylogeny/allB_veryfasttree.nwk"
    threads: 32
    log: f"{LOGS}/06_veryfasttree.log"
    conda: "../envs/veryfasttree.yaml"
    shell:
        "VeryFastTree -gtr -nt -threads {threads} {input.var} > {output.nwk} 2> {log}"
