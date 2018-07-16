get_outdir = lambda w, output: os.path.dirname(output[0])


rule ustacks:
    input:
        "trimmed/{individual}.R1.fastq.gz",
        "trimmed/{individual}.R2.fastq.gz"
    output:
        "ustacks/{individual}.M={max_individual_mm}.m={min_reads}/{individual}.tags.tsv",
        "ustacks/{individual}.M={max_individual_mm}.m={min_reads}/{individual}.snps.tsv",
        "ustacks/{individual}.M={max_individual_mm}.m={min_reads}/{individual}.alleles.tsv"
    params:
        outdir=get_outdir
    threads: 8
    conda:
        "../envs/stacks.yaml"
    shell:
        "ustacks -p {threads} -f {input} -o {output} "
        "--name {wildcards.individual} "
        "-M {wildcards.max_individual_mm} "
        "-m {wildcards.min_reads}"


def fmt_ustacks_input(wildcards, input):
    return ["-s {}".format(os.path.dirname(f)) for f in input.ustacks]

ustacks_individuals = expand(
    "ustacks/{individual}.M={{max_individual_mm}}.m={{min_reads}}/{individual}.tags.tsv",
    individual=individuals.id)


rule cstacks:
    input:
        ustacks=ustacks_individuals
    output:
        "cstacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/tags.tsv",
        "cstacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/snps.tsv",
        "cstacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/alleles.tsv"
    params:
        outdir=get_outdir,
        individuals=fmt_ustacks_input
    conda:
        "../envs/stacks.yaml"
    threads: 8
    shell:
        "cstacks -p {threads} {params.individuals} -o {params.outdir}"


rule sstacks:
    input:
        ustacks=ustacks_individuals,
        cstacks=rules.cstacks.output[0]
    output:
        expand("sstacks/n={{max_locus_mm}}.M={{max_individual_mm}}.m={{min_reads}}/{individual}.matches.tsv",
               individual=individuals.id)
    params:
        outdir=get_outdir,
        individuals=fmt_ustacks_input
    conda:
        "../envs/stacks.yaml"
    threads: 8
    shell:
        "sstacks -p {threads} {params.individuals} -c {input.cstacks} "
        "-o {params.outdir}"


rule tsv2bam:
    input:
        sstacks=rules.sstacks.output,
        reads=["trimmed/{individual}.R1.fastq.gz",
               "trimmed/{individual}.R2.fastq.gz"]
    output:
        "sstacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/{individual}.bam"
    params:
        sstacks_dir=lambda w, output: os.path.dirname(output[0]),
        read_dir=lambda w, input: os.path.dirname(input.reads[0])
    conda:
        "../envs/stacks.yaml"
    shell:
        "tsv2bam -s {wildcards.individual} -R {params.read_dir} "
        "-P {params.sstacks_dir}"


rule gstacks:
    input:
        expand("sstacks/n={{max_locus_mm}}.M={{max_individual_mm}}.m={{min_reads}}/{individual}.bam",
               individual=individuals.id)
    output:
        "gstacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/calls.vcf"
    params:
        outdir=get_outdir,
        bams=lambda w, input: ["-B {}".format(f) for f in input],
        config=config["params"]["gstacks"]
    conda:
        "../envs/stacks.yaml"
    threads: 8
    shell:
        "gstacks {params.config} {params.bams} -O {params.outdir}"
