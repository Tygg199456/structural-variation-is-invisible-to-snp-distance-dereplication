# -------------------------------------------------------------- 03_ancient ---
# Tier C. Storage is the binding constraint, so this runs in two phases:
#   phase 1 (pilot): depth >= 10x, ~30 libraries, 200-400 GB
#   phase 2 (full) : depth >=  5x, ~60-80 libraries, 0.5-2 TB
# Run the pilot end to end before committing the disk for phase 2.
#
# Tier C may contribute SNPs, deep calibration points, and COARSE presence/
# absence (breadth >= 90% at >= 3x). It may NEVER contribute frameshift or
# pseudogene calls: at these coverages an indel call is indistinguishable
# from a damage artefact or a stochastic gap.

rule amdirt_export:
    """
    AncientMetagenomeDir / AMDirT -> ENA run accessions + library metadata.
    Library prep matters: UDG-treated libraries carry little or no terminal
    deamination and will fail the damage gate. Record prep and exempt them.
    """
    output:
        samples = f"{DATA}/ancient/amdirt_samples.tsv",
        libs    = f"{DATA}/ancient/amdirt_libraries.tsv"
    log: f"{LOGS}/03_amdirt.log"
    conda: "../envs/amdirt.yaml"
    shell:
        r"""
        AMDirT filter --table ancientmetagenome-hostassociated \
          --query "species == 'Yersinia pestis'" \
          --output $(dirname {output.samples}) &> {log}
        """


rule eager_run:
    """
    nf-core/eager: adapter removal, aDNA-tuned bwa aln, mapDamage2
    authentication, genotyping. Standardised and reviewer-legible.
    """
    input:
        sheet = f"{DATA}/ancient/eager_samplesheet.tsv",
        ref   = f"{DATA}/ref/{REF}.fna"
    output:
        done = f"{RES}/ancient/eager/.done"
    threads: 32
    log: f"{LOGS}/03_eager.log"
    conda: "../envs/nextflow.yaml"
    shell:
        r"""
        nextflow run nf-core/eager -r 2.5.0 -profile conda \
          --input {input.sheet} --fasta {input.ref} \
          --mapper bwaaln --bwaalnn 0.01 --bwaalno 2 \
          --run_bam_filtering --bam_mapping_quality_threshold 30 \
          --run_damage_calculation \
          --outdir {RES}/ancient/eager &> {log}
        touch {output.done}
        """


checkpoint tier_c_select:
    input:
        done = f"{RES}/ancient/eager/.done",
        libs = f"{DATA}/ancient/amdirt_libraries.tsv"
    output:
        keep = f"{RES}/ancient/tierC_keep.tsv"
    params:
        min_depth = config["tiers"]["C"]["pilot_min_depth"],
        auth      = config["tiers"]["C"]["authentication"]
    log: f"{LOGS}/03_tierC_select.log"
    conda: "../envs/py.yaml"
    script: "../scripts/tier_c_select.py"


rule coarse_pav_ancient:
    """Breadth-based gene presence for ancient genomes. No frameshift calls."""
    input:
        bams = lambda w: expand(f"{RES}/ancient/eager/deduplication/{{s}}/{{s}}.bam",
                                s=tier_c()),
        gff  = f"{RES}/annotate/{REF}/{REF}.gff3"
    output:
        tsv = f"{RES}/ancient/coarse_pav.tsv"
    params:
        min_breadth = config["tiers"]["C"]["coarse_pav"]["min_breadth"],
        min_depth   = config["tiers"]["C"]["coarse_pav"]["min_depth"]
    log: f"{LOGS}/03_coarse_pav.log"
    conda: "../envs/py.yaml"
    script: "../scripts/coarse_pav.py"


rule figS2_downsample_control:
    """
    Fig.S2 -- the control a reviewer will demand.
    Downsample high-coverage Tier A genomes to Tier C depth, rerun the same
    presence/absence logic, and report the false-loss rate. Without this the
    ancient gene losses in Fig.3 are not defensible.
    """
    input:
        tierA   = f"{RES}/dedup/tierA_keep.tsv",
        depths  = f"{RES}/ancient/tierC_keep.tsv"
    output:
        pdf = f"{RES}/figures/figS2_downsample_control.pdf",
        tsv = f"{RES}/tables/figS2_false_loss_rate.tsv"
    threads: 16
    log: f"{LOGS}/03_downsample_control.log"
    conda: "../envs/py.yaml"
    script: "../scripts/downsample_control.py"
