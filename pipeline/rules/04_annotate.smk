# ------------------------------------------------------------- 04_annotate ---
# One Bakta DB version for every genome, references included.
# Mixing PGAP vintages is the largest single source of spurious pseudogene
# differences between genomes, and it is invisible downstream.

rule bakta:
    input:
        fa = f"{DATA}/genomes/{{acc}}/{{acc}}.fna"
    output:
        gff = f"{RES}/annotate/{{acc}}/{{acc}}.gff3",
        faa = f"{RES}/annotate/{{acc}}/{{acc}}.faa",
        tsv = f"{RES}/annotate/{{acc}}/{{acc}}.tsv"
    threads: 8
    params:
        db = config["annotation"]["db_version"]
    log: f"{LOGS}/04_bakta/{{acc}}.log"
    conda: "../envs/bakta.yaml"
    shell:
        r"""
        bakta --db $BAKTA_DB --threads {threads} --force \
          --prefix {wildcards.acc} \
          --output $(dirname {output.gff}) \
          --genus Yersinia --species pestis \
          --keep-contig-headers {input.fa} &> {log}
        """


rule pseudofinder:
    """
    Tier A only. Reported calls require agreement between pseudofinder and
    bakta's own `pseudo` tag, which trades recall for a far lower false
    positive rate -- the right side of that trade when the pseudogene count
    IS the measurement.
    """
    input:
        gbk = f"{RES}/annotate/{{acc}}/{{acc}}.gbff",
        db  = f"{DATA}/db/yersinia_enterobacterales.dmnd"
    output:
        tsv = f"{RES}/pseudogene/{{acc}}/{{acc}}_pseudos.gff"
    threads: 8
    log: f"{LOGS}/04_pseudofinder/{{acc}}.log"
    conda: "../envs/pseudofinder.yaml"
    shell:
        r"""
        pseudofinder.py annotate --genome {input.gbk} \
          --database {input.db} --diamond --threads {threads} \
          --outprefix {RES}/pseudogene/{wildcards.acc}/{wildcards.acc} &> {log}
        """


rule pseudo_matrix:
    """Tier A pseudo calls, agreement-filtered.

    Emits a LONG table, not a wide matrix. Coordinates are strain-specific, so
    a strain x coordinate matrix would have columns that are not comparable
    across rows -- a file shaped like a matrix eventually gets used as one.
    The real strain x ortholog matrix is built by combine_pav.py on Panaroo
    orthogroups.
    """
    input:
        gffs  = lambda w: expand(f"{RES}/pseudogene/{{acc}}/{{acc}}_pseudos.gff",
                                 acc=tier_a()),
        bakta = lambda w: expand(f"{RES}/annotate/{{acc}}/{{acc}}.tsv",
                                 acc=tier_a())
    output:
        long   = f"{RES}/pseudogene/pseudo_calls_long.tsv",
        counts = f"{RES}/pseudogene/pseudo_counts.tsv"
    params:
        require_agreement = config["pseudogene"]["require_bakta_agreement"]
    log: f"{LOGS}/04_pseudo_matrix.log"
    conda: "../envs/py.yaml"
    script: "../scripts/pseudo_matrix.py"


rule isescan:
    input:
        fa = f"{DATA}/genomes/{{acc}}/{{acc}}.fna"
    output:
        tsv = f"{RES}/is_elements/{{acc}}/{{acc}}.fna.tsv"
    threads: 4
    log: f"{LOGS}/04_isescan/{{acc}}.log"
    conda: "../envs/isescan.yaml"
    shell:
        "isescan.py --seqfile {input.fa} --output {RES}/is_elements/{wildcards.acc} "
        "--nthread {threads} &> {log}"


rule is_breakpoint_overlap:
    """IS copies vs loss/pseudogenisation breakpoints -- the mechanism panel."""
    input:
        is_tsv = lambda w: expand(f"{RES}/is_elements/{{acc}}/{{acc}}.fna.tsv",
                                  acc=tier_a()),
        pav    = f"{RES}/pangenome/gene_presence_absence.Rtab"
    output:
        tsv = f"{RES}/is_elements/breakpoint_overlap.tsv"
    log: f"{LOGS}/04_is_breakpoints.log"
    conda: "../envs/py.yaml"
    script: "../scripts/is_breakpoints.py"
