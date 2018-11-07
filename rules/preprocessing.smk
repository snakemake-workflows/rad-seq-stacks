
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


rule mark_duplicates:
    input:
        fq1=lambda w: units.loc[w.unit, "fq1"],
        fq2="trimmed-spacer/{unit}.2.fq.gz"
    output:
        "dedup/{unit}.markdup.bam"
    log:
        "logs/mark-duplicates/{unit}.log"
    conda:
        "../envs/markdup.yaml"
    params:
        dbr_len=config["dbr"]["len"],
        dbr_dist=config["dbr"]["max_dist"],
        seq_dist=config["dbr"]["max_seq_dist"]
    script:
        "../scripts/mark-duplicates.py"


rule generate_consensus_reads:
    input:
        "dedup/{unit}.markdup.bam"
    output:
        consensus="dedup/{unit}.consensus.bam",
        singletons="dedup/{unit}.singletons.bam",
    conda:
        "../envs/fgbio.yaml"
    shell:
        "fgbio CallMolecularConsensusReads --min-input-base-quality 0 "
        "--input {input} --output {output.consensus} "
        "--min-reads 2 --rejects {output.singletons}"


rule bam_to_fastq:
    input:
        "dedup/{unit}.consensus.bam",
        "dedup/{unit}.singletons.bam"
    output:
        "dedup/{unit}.consensus.1.fq.gz",
        "dedup/{unit}.consensus.2.fq.gz"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        for f in {input}
        do
            samtools fastq -1 {output[0]} -2 {output[1]} $f
        done
        """


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
            "-a {}".format(config["adapter"]) if config.get("adapter") else "")
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
