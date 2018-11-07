
rule barcodes:
    output:
        "barcodes/{unit}.tsv"
    run:

        d = individuals.loc[individuals.unit == wildcards.unit, ["p5_barcode", "id"]]
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


ruleorder: trim_p7_spacer > zip_fastq


rule unzip_fastq:
    input:
        "{prefix}.{ext}.gz"
    output:
        temp("{prefix}.{ext,fastq|fq}")
    shell:
        "gzip -d -c {input} > {output}"


rule zip_fastq:
    input:
        "{prefix}.fastq"
    output:
        "{prefix}.fq.gz"
    shell:
        "gzip -c {input} > {output}"


def calib_fq_input(wildcards):
    fq1 = units.loc[wildcards.unit, "fq1"]
    if fq1.endswith(".gz"):
        fq1 = fq1[:-3]
    return fq1, "trimmed-spacer/{unit}.2.fq".format(**wildcards)


rule calib_cluster:
    input:
        calib_fq_input
    output:
        "dedup/{unit}.cluster"
    log:
        "logs/calib/{unit}.log"
    params:
        dbr_len=config["dbr"]["len"],
        prefix=lambda w, output: output[0][:-7]
    shell:
        "calib -f {input[0]} -r {input[1]} "
        "-l {params.dbr_len} -o {params.prefix} > {log}"


rule group_by_dbr_cluster:
    input:
        "dedup/{unit}.cluster"
    output:
        "dedup/{unit}.dbr-grouped.bam"
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/group-by-dbr-cluster.py"


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


# rule calib_consensus:
#     input:
#         fq=calib_fq_input,
#         cluster="dedup/{unit}.cluster"
#     output:
#         temp("dedup/{unit}.consensus.1.fastq"),
#         temp("dedup/{unit}.consensus.2.fastq")
#     params:
#         prefix=lambda w, output: output[0][:-8]
#     shell:
#         "calib_cons -c {input.cluster} -q {input.fq[0]} {input.fq[1]} "
#         "-o {params.prefix}.1 {params.prefix}.2"


rule extract:
    input:
        fq1=expand("dedup/{unit}.consensus.1.fq.gz", unit=units.id),
        fq2=expand("dedup/{unit}.consensus.2.fq.gz", unit=units.id),
        barcodes=expand("barcodes/{unit}.tsv", unit=units.id)
    output:
        expand(["extracted/{individual}.1.fq.gz",
                "extracted/{individual}.2.fq.gz"],
               individual=individuals.id)
    log:
        expand("logs/extract/{individual}.log",
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
