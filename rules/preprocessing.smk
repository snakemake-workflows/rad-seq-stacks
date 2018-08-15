
rule barcodes:
    output:
        "barcodes/{unit}.tsv"
    run:
        d = individuals[["p5_barcode", "id"]]
        #d["p7_barcode"] = units.loc[wildcards.unit, "p7_barcode"]
        d[["p5_barcode", "id"]].to_csv(output[0], index=False, header=None, sep="\t")


rule trim_spacer:
    input:
        lambda w: units.loc[w.unit, "fq2"]
    output:
        "trimmed-spacer/{unit}.2.fq.gz"
    conda:
        "../envs/seqtk.yaml"
    shell:
        # TODO look up spacer length in units.tsv
        "seqtk trimfq -b3 {input} > {output}"


rule extract:
    input:
        fq1=units.fq1,
        fq2=expand("trimmed-spacer/{unit}.2.fq.gz", unit=units.id),
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


rule remove_duplicates:
    input:
        "extracted/{individual}.1.fq.gz",
        "extracted/{individual}.2.fq.gz"
    output:
        "nodup/{individual}.1.fq.gz",
        "nodup/{individual}.2.fq.gz"
    shell:
        # TODO remove trimming, properly handle DBR to remove duplicates
        "cp {input[0]} {output[0]}; "
        "cp {input[1]} {output[1]};"


rule trim:
    input:
        "nodup/{individual}.1.fq.gz",
        "nodup/{individual}.2.fq.gz"
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
