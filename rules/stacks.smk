rule ustacks:
    input:
        "trimmed/{individual}.1.fq.gz"
    output:
        "ustacks/M={max_individual_mm}.m={min_reads}/{individual}.tags.tsv.gz",
        "ustacks/M={max_individual_mm}.m={min_reads}/{individual}.snps.tsv.gz",
        "ustacks/M={max_individual_mm}.m={min_reads}/{individual}.alleles.tsv.gz"
    params:
        outdir=get_outdir,
        hash=lambda w: individuals.loc[w.individual, "hash"]
    threads: 8
    conda:
        "../envs/stacks.yaml"
    shell:
        "ustacks -p {threads} -f {input} -o {params.outdir} "
        "--name {wildcards.individual} "
        "-i {params.hash} "
        "-M {wildcards.max_individual_mm} "
        "-m {wildcards.min_reads}"


def fmt_ustacks_input(wildcards, input):
    return ["-s {}".format(f[:-len(".tags.tsv.gz")]) for f in input.ustacks]

ustacks_individuals = expand(
    "ustacks/M={{max_individual_mm}}.m={{min_reads}}/{individual}.tags.tsv.gz",
    individual=individuals.id)


rule cstacks:
    input:
        ustacks=ustacks_individuals
    output:
        "stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.tags.tsv.gz",
        "stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.snps.tsv.gz",
        "stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.alleles.tsv.gz"
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
        expand("stacks/n={{max_locus_mm}}.M={{max_individual_mm}}.m={{min_reads}}/{individual}.matches.tsv.gz",
               individual=individuals.id),
    params:
        outdir=get_outdir,
        individuals=fmt_ustacks_input,
        cstacks_dir=lambda w, input: os.path.dirname(input.cstacks)
    conda:
        "../envs/stacks.yaml"
    threads: 8
    shell:
        "sstacks -p {threads} {params.individuals} -c {params.cstacks_dir} "
        "-o {params.outdir}"


rule link_ustacks:
    input:
        "ustacks/M={max_individual_mm}.m={min_reads}/{individual}.{type}.tsv.gz"
    output:
        "stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/{individual}.{type}.tsv.gz",
    shell:
        "ln -s -r {input} {output}"


rule tsv2bam:
    input:
        sstacks=rules.sstacks.output,
        ustacks=expand("stacks/n={{max_locus_mm}}.M={{max_individual_mm}}.m={{min_reads}}/{{individual}}.{type}.tsv.gz",
                       type=["tags", "snps", "alleles"]),
        reads=["trimmed/{individual}.1.fq.gz",
               "trimmed/{individual}.2.fq.gz"]
    output:
        "stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/{individual}.matches.bam"
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
        bams=expand("stacks/n={{max_locus_mm}}.M={{max_individual_mm}}.m={{min_reads}}/{individual}.matches.bam",
                    individual=individuals.id),
        popmap="population-map.tsv"
    output:
        "stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.calls",
        "stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.fa.gz"
    params:
        outdir=get_outdir,
        bam_dir=lambda w, input: os.path.dirname(input.bams[0]),
        config=config["params"]["gstacks"]
    conda:
        "../envs/stacks.yaml"
    threads: 8
    shell:
        "gstacks {params.config} -P {params.bam_dir} -O {params.outdir} "
        "-M {input.popmap}"


rule populations:
    input:
        "stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.calls"
    output:
        "calls/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/populations.snps.vcf"
    params:
        outdir=get_outdir,
        gstacks_dir=lambda w, input: os.path.dirname(input[0])
    conda:
        "../envs/stacks.yaml"
    threads: 8
    shell:
        "populations -t {threads} -P {params.gstacks_dir} -O {params.outdir} --vcf"
