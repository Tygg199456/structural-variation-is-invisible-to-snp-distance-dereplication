# ---------------------------------------------------------------- 00_fetch ---
# Pull assembly summaries + sequences from NCBI Datasets, and BioSample metadata.

rule datasets_summary_all:
    output:
        jsonl = f"{DATA}/ncbi/all_assemblies.jsonl"
    params:
        taxid = config["project"]["taxid"]
    log: f"{LOGS}/00_datasets_summary_all.log"
    conda: "../envs/ncbi.yaml"
    shell:
        r"""
        datasets summary genome taxon {params.taxid} --as-json-lines \
          > {output.jsonl} 2> {log}
        wc -l {output.jsonl} >> {log}
        """


rule datasets_summary_complete:
    output:
        jsonl = f"{DATA}/ncbi/complete_chromosome.jsonl"
    params:
        taxid = config["project"]["taxid"]
    log: f"{LOGS}/00_datasets_summary_complete.log"
    conda: "../envs/ncbi.yaml"
    shell:
        r"""
        datasets summary genome taxon {params.taxid} \
          --assembly-level complete,chromosome --as-json-lines \
          > {output.jsonl} 2> {log}
        """


rule flatten_metadata:
    """
    JSONL -> flat TSV, with missing-value tokens normalised.
    Reproduces the measured completeness numbers recorded in config.yaml.
    """
    input:
        allj  = f"{DATA}/ncbi/all_assemblies.jsonl",
        compj = f"{DATA}/ncbi/complete_chromosome.jsonl"
    output:
        tsv  = f"{RES}/qc/assembly_metadata.tsv",
        comp = f"{RES}/qc/metadata_completeness.tsv"
    params:
        missing = config["metadata"]["missing_tokens"]
    log: f"{LOGS}/00_flatten_metadata.log"
    conda: "../envs/py.yaml"
    script: "../scripts/flatten_metadata.py"


rule download_assemblies:
    """
    Fetch FASTA for one accession. Downloads are untrusted input: each lands
    in its own directory and is never executed or imported from.
    """
    output:
        fa = f"{DATA}/genomes/{{acc}}/{{acc}}.fna"
    log: f"{LOGS}/00_download/{{acc}}.log"
    retries: 3
    conda: "../envs/ncbi.yaml"
    shell:
        r"""
        tmp=$(mktemp -d)
        datasets download genome accession {wildcards.acc} \
          --include genome --filename $tmp/{wildcards.acc}.zip &> {log}
        unzip -q -o $tmp/{wildcards.acc}.zip -d $tmp/x
        mkdir -p $(dirname {output.fa})
        cat $tmp/x/ncbi_dataset/data/{wildcards.acc}/*.fna > {output.fa}
        rm -rf $tmp
        """


rule fetch_reference:
    output:
        fa  = f"{DATA}/ref/{REF}.fna",
        gbk = f"{DATA}/ref/{REF}.gbff"
    log: f"{LOGS}/00_fetch_reference.log"
    conda: "../envs/ncbi.yaml"
    shell:
        r"""
        tmp=$(mktemp -d)
        datasets download genome accession {REF} \
          --include genome,gbff --filename $tmp/ref.zip &> {log}
        unzip -q -o $tmp/ref.zip -d $tmp/x
        cat $tmp/x/ncbi_dataset/data/{REF}/*.fna  > {output.fa}
        cat $tmp/x/ncbi_dataset/data/{REF}/*.gbff > {output.gbk}
        rm -rf $tmp
        """
