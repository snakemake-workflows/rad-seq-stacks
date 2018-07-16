get_outdir = lambda w, output: os.path.dirname(output[0])


rule ustacks:
    input:
        "trimmed/{individual}.R1.fastq.gz",
        "trimmed/{individual}.R2.fastq.gz"
    output:
        "ustacks/{individual}.M={max_individual_mm}.m={max_reads}/{individual}.tags.tsv",
        "ustacks/{individual}.M={max_individual_mm}.m={max_reads}/{individual}.snps.tsv",
        "ustacks/{individual}.M={max_individual_mm}.m={max_reads}/{individual}.alleles.tsv"
    params:
        outdir=get_outdir
    threads: 8
    shell:
        "ustacks -p {threads} -f {input} -o {output} "
        "--name {wildcards.individual} "
        "-M {wildcards.max_mismatches} "
        "-m {wildcards.max_reads}"


def fmt_ustacks_input(wildcards, input):
    return ["-s {}".format(os.path.dirname(f)) for f in input.ustacks]


rule cstacks:
    input:
        ustacks=expand(rules.ustacks.output[0], individual=individuals.id)
    output:
        "cstacks/n={max_locus_mm}.M={max_individual_mm}.m={max_reads}/tags.tsv",
        "cstacks/n={max_locus_mm}.M={max_individual_mm}.m={max_reads}/snps.tsv",
        "cstacks/n={max_locus_mm}.M={max_individual_mm}.m={max_reads}/alleles.tsv"
    params:
        outdir=get_outdir,
        individuals=fmt_ustacks_input
    threads: 8
    shell:
        "cstacks -p {threads} {params.individuals} -o {params.outdir}"


rule sstacks:
    input:
        ustacks=expand(rules.ustacks.output[0], individual=individuals.id),
        cstacks=rules.cstacks.output[0]
    output:
        "sstacks/n={max_locus_mm}.M={max_individual_mm}.m={max_reads}/{individuals}.matches.tsv"
    params:
        outdir=get_outdir,
        individuals=fmt_ustacks_input
    threads: 8
    shell:
        "sstacks -p {threads} {params.individuals} -c {input.cstacks} "
        "-o {params.outdir}"


rule tsv2bam:
    input:
        sstacks=rules.sstacks.output[0],
        reads=["trimmed/{individual}.R1.fastq.gz",
               "trimmed/{individual}.R2.fastq.gz"]
    output:
        "sstacks/n={max_locus_mm}.M={max_individual_mm}.m={max_reads}/{individual}.bam"
    params:
        sstacks_dir=lambda w, output: os.path.dirname(output),
        read_dir=lambda w, input: os.path.dirname(input.reads[0])
    shell:
        "tsv2bam -s {wildcard.individual} -R {params.read_dir} "
        "-P {params.sstacks}"


rule gstacks:
    input:
        rules.tsv2bam.output
    output:
        "gstacks/n={max_locus_mm}.M={max_individual_mm}.m={max_reads}/calls.vcf"
    params:
        outdir=get_outdir,
        bams=lambda w, input: ["-B {}".format(f) for f in input],
        config=config["params"]["gstacks"]
    threads: 8
    shell:
        "gstacks {params.config} {params.bams} -O {params.outdir}"
