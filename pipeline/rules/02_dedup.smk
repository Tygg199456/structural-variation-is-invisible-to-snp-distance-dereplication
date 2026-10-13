# ---------------------------------------------------------------- 02_dedup ---
# Deduplication by core-SNP distance, computed from per-sample VCFs.
#
# Two things this deliberately does NOT do:
#   1. ANI. Every Y. pestis pair exceeds 99.99% ANI, so ANI cannot separate
#      strains in this species. NCBI's best_ani_match.ani compares an assembly
#      to the species type strain, not to the other assemblies.
#   2. Build a whole-genome alignment first. 3453 x 4.65 Mb would be ~16 GB.
#      The genotype matrix below is 3453 x ~25k int8 = ~86 MB.

rule snippy_single:
    input:
        asm = f"{DATA}/genomes/{{acc}}/{{acc}}.fna",
        ref = f"{DATA}/ref/{REF}.gbff"
    output:
        vcf = f"{RES}/snippy/{{acc}}/snps.vcf",
        aln = f"{RES}/snippy/{{acc}}/snps.aligned.fa"
    threads: 4
    params:
        mincov  = config["variants"]["min_cov"],
        minfrac = config["variants"]["min_frac"]
    log: f"{LOGS}/02_snippy/{{acc}}.log"
    conda: "../envs/snippy.yaml"
    shell:
        r"""
        snippy --force --cpus {threads} \
          --outdir $(dirname {output.vcf}) \
          --ref {input.ref} --ctgs {input.asm} \
          --mincov {params.mincov} --minfrac {params.minfrac} &> {log}
        """


rule genotype_matrix:
    """Union of variant sites across one set -> int8 matrix (~86 MB Tier B).

    Parameterised by {set} so Tier A (256) and Tier B (3,441) build separate
    matrices. Without this, the Tier A dedup would transitively require a
    snippy run for every Tier B accession, making Gate 1c cost Tier B money.
    snipsy output is per-sample and cached, so Tier A's runs are reused when
    Tier B is built.
    """
    input:
        vcfs = lambda w: expand(f"{RES}/snippy/{{acc}}/snps.vcf",
                                acc=(tier_a_candidates() if w.set == "tierA"
                                     else tier_b_candidates()))
    output:
        npz   = f"{RES}/dedup/genotypes_{{set}}.npz",
        sites = f"{RES}/dedup/sites_{{set}}.tsv"
    threads: 8
    log: f"{LOGS}/02_genotype_matrix_{{set}}.log"
    conda: "../envs/py.yaml"
    script: "../scripts/genotype_matrix.py"


rule snp_distance:
    """
    All-pairs SNP distance via one-hot encoding + BLAS matmul.
    ~3453 x 3453 in roughly 30 s; no pairwise Python loop.

    Emits BOTH matrices. The raw one must be a declared output: v5 wrote it as
    an undeclared side effect, so snakemake saw the downstream consumer's input
    as unproducible and refused to build the DAG (gate1d died in the first
    second with MissingInputException). The threshold applies to raw.
    """
    input:
        npz = f"{RES}/dedup/genotypes_{{set}}.npz"
    output:
        dist = f"{RES}/dedup/snp_distance_{{set}}.npy",
        raw  = f"{RES}/dedup/snp_distance_raw_{{set}}.npy"
    threads: 16
    log: f"{LOGS}/02_snp_distance_{{set}}.log"
    conda: "../envs/py.yaml"
    script: "../scripts/snp_distance.py"


checkpoint tier_a_dedup:
    """Tier A: complete-linkage clustering at the calibrated threshold.

    Reads the RAW matrix, not the coverage-scaled one: merging must rest on
    observed differences. At t=0 the two are identical, which is why using the
    scaled matrix went unnoticed until same-strain pairs were measured 3-21
    apart and A1122 (raw 18, scaled 21) failed to merge at t=20.
    """
    input:
        dist    = f"{RES}/dedup/snp_distance_raw_tierA.npy",
        cand    = f"{RES}/qc/tierA_refined.tsv",
        meta    = f"{RES}/qc/assembly_metadata.tsv",
        profile = f"{RES}/qc/deletion_profile_all.tsv"
    output:
        keep     = f"{RES}/dedup/tierA_keep.tsv",
        clusters = f"{RES}/dedup/tierA_clusters.tsv"
    params:
        max_dist  = config["dedup"]["max_snp_distance"],
        priority  = config["dedup"]["representative_priority"],
        whitelist = config["whitelist"]
    log: f"{LOGS}/02_tierA_dedup.log"
    conda: "../envs/py.yaml"
    script: "../scripts/snp_dedup.py"


checkpoint tier_b_backbone:
    """
    Tier B backbone: same dedup, plus a hard cap on taxa.
    Rationale for capping rather than using all 3453:
      - 678 genomes already yield only ~23k non-repetitive core SNPs over
        3.8 Mb, so most terminal branches carry single-digit SNP counts and
        many internal branches are zero-length.
      - UFBoot support on zero-length branches is noise, not information.
      - Madagascar and Brazil are heavily over-sampled in public archives;
        leaving them unthinned lets them dictate the topology.
      - Strict-core site retention falls monotonically with each gappy draft
        added, so more taxa means fewer usable sites.
    Collapsed members are retained as cluster sizes and drawn as bubbles on
    the Fig.1 outer ring.
    """
    input:
        dist    = f"{RES}/dedup/snp_distance_raw_tierB.npy",
        cand    = f"{RES}/qc/tierB_candidates.tsv",
        meta    = f"{RES}/qc/assembly_metadata.tsv",
        profile = f"{RES}/qc/deletion_profile_all.tsv"
    output:
        keep     = f"{RES}/dedup/tierB_backbone.tsv",
        clusters = f"{RES}/dedup/tierB_clusters.tsv"
    params:
        max_dist  = config["dedup"]["max_snp_distance"],
        hard_cap  = config["tiers"]["B"]["backbone_hard_cap"],
        priority  = config["dedup"]["representative_priority"],
        whitelist = config["whitelist"]
    log: f"{LOGS}/02_tierB_backbone.log"
    conda: "../envs/py.yaml"
    script: "../scripts/snp_dedup.py"
