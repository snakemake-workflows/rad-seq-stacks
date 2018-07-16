def get_fastq(wildcards):
    return units.loc[individuals.loc[wildcards.individual, "unit"], ["fq1", "fq2"]]


rule barcodes:
    output:
        "barcodes/{unit}.tsv"
    run:
        individuals[["barcode", "id"]].to_csv(output[0], header=None)



rule extract:
    input:
        reads=get_fastq,
        barcodes="barcodes/{unit}.tsv"
    output:
        "extracted/{individual}.R1.fastq.gz",
        "extracted/{individual}.R2.fastq.gz"
    params:
        enzymes=config["restriction-enzyme"]
    shell:
        "process_radtags -1 {input[0]} -2 {input[0]} "
        "--renz_1 {params.enzymes[0] --renz_2 {params.enzymes[1]} "
        "-b {input.barcodes}"


rule remove_duplicates:
    input:
        "extracted/{individual}.R1.fastq.gz",
        "extracted/{individual}.R2.fastq.gz"
    output:
        "nodup/{individual}.R1.fastq.gz",
        "nodup/{individual}.R2.fastq.gz"
    shell:
        "cp {input[0]} {output[0]}; cp {input[1]} {output[1]}; "


rule trim:
    input:
        "nodup/{individual}.R1.fastq.gz",
        "nodup/{individual}.R2.fastq.gz"
    output:
        fastq1="trimmed/{individual}.R1.fastq.gz",
        fastq2="trimmed/{individual}.R2.fastq.gz",
        qc="trimmed/{individual}.qc.txt"
    wrapper:
        "0.27.1/bio/cutadapt/pe"
