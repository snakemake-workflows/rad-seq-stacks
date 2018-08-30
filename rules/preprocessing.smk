
rule barcodes:
    output:
        "barcodes/{unit}.tsv"
    run:
        d = individuals[["p5_barcode", "id"]]
        #d["p7_barcode"] = units.loc[wildcards.unit, "p7_barcode"]
        d[["p5_barcode", "id"]].to_csv(output[0], index=False, header=None, sep="\t")


rule trim_p7_spacer:
    input:
        lambda w: units.loc[w.unit, "fq2"]
    output:
        "trimmed-spacer/{unit}.2.fq.gz"
    params:
        spacer=lambda w: units.loc[w.unit, "p7_spacer"]
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk trimfq -b {params.spacer} {input} | gzip > {output}"


rule fastq_to_bam:
    input:
        lambda w: units.loc[w.unit, "fq1"],
        "trimmed-spacer/{unit}.2.fq.gz"
    output:
        "dedup/{unit}.dbr-annotated.bam"
    params:
        dbr=config["dbr"]["len"]
    conda:
        "../envs/fgbio.yaml"
    shell:
        "fgbio FastqToBam --input {input} --output {output} "
        "--sample {wildcards.unit} --library rad-seq "
        "--read-structures +T {params.dbr}M+T"


rule fake_mapping:
    input:
        "dedup/{unit}.dbr-annotated.bam"
    output:
        "dedup/{unit}.fake-mapping.bam"
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/fake-mapping.py"


rule group_by_dbr:
    input:
        "dedup/{unit}.fake-mapping.bam"
    output:
        "dedup/{unit}.dbr-grouped.bam"
    params:
        strategy=config["dbr"]["grouping-strategy"]
    conda:
        "../envs/fgbio.yaml"
    shell:
        #"umi_tools group -I {input} --paired --output-bam -S {output} "
        #"--umi-group-tag MI"
        "fgbio GroupReadsByUmi --input {input} --output {output} "
        "--strategy {params.strategy} --min-map-q 0 --include-non-pf-reads"



rule generate_consensus_reads:
    input:
        "dedup/{unit}.dbr-grouped.bam"
    output:
        "dedup/{unit}.consensus.bam"
    conda:
        "../envs/fgbio.yaml"
    shell:
        "fgbio CallMolecularConsensusReads --input {input} --output {output} "
        "--min-reads 1"


rule bam_to_fastq:
    input:
        "dedup/{unit}.consensus.bam"
    output:
        "dedup/{unit}.consensus.1.fq.gz",
        "dedup/{unit}.consensus.2.fq.gz"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools fastq -1 {output[0]} -2 {output[1]} {input}"


rule extract:
    input:
        fq1=expand("dedup/{unit}.consensus.1.fq.gz", unit=units.id),
        fq2=expand("dedup/{unit}.consensus.2.fq.gz", unit=units.id),
        barcodes=expand("barcodes/{unit}.tsv", unit=units.id)
    output:
        expand(["extracted/{individual}.1.fq.gz",
                "extracted/{individual}.2.fq.gz"],
               individual=individuals.id)
    params:
        enzymes=config["restriction-enzyme"],
        outdir=get_outdir,
        units=units,
        extra=config["params"]["process_radtags"]
    conda:
        "../envs/stacks.yaml"
    script:
        "../scripts/extract-individuals.py"


rule trim:
    input:
        "extracted/{individual}.1.fq.gz",
        "extracted/{individual}.2.fq.gz"
    output:
        fastq1="trimmed-adapter/{individual}.1.fq.gz",
        fastq2="trimmed-adapter/{individual}.2.fq.gz",
        qc="trimmed/{individual}.qc.txt"
    params:
        config["params"]["cutadapt"] + (
            "-a {}".format(config["adapter"]) if config["adapter"] else "")
    wrapper:
        "0.27.1/bio/cutadapt/pe"


rule force_same_length:
    input:
        "trimmed-adapter/{individual}.{read}.fq.gz"
    output:
        "trimmed/{individual}.{read}.fq.gz"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "len=`seqtk fqchk {input} | grep -oP 'min_len: \\K[0-9]+'`; "
        "seqtk trimfq -L$len {input} | gzip -c > {output}"


rule population_map:
    output:
        "population-map.tsv"
    run:
        d = individuals[["id"]]
        d["pop"] = 1
        d.to_csv(output[0], index=False, header=None, sep="\t")
