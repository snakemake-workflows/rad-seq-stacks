
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
        "seqtk trimfq -L91 {input[0]} | gzip -c > {output[0]}; "
        "seqtk trimfq -L91 {input[1]} | gzip -c > {output[1]}; "


rule trim:
    input:
        "nodup/{individual}.1.fq.gz",
        "nodup/{individual}.2.fq.gz"
    output:
        fastq1="trimmed/{individual}.1.fq.gz",
        fastq2="trimmed/{individual}.2.fq.gz",
        qc="trimmed/{individual}.qc.txt"
    wrapper:
        "0.27.1/bio/cutadapt/pe"


rule population_map:
    output:
        "population-map.tsv"
    run:
        d = individuals[["id"]]
        d["pop"] = 1
        d.to_csv(output[0], index=False, header=None, sep="\t")
