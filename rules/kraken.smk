rule kraken:
    input:
        reads=lambda w: units.loc[w.unit, ["fq1", "fq2"]],
        db=config["params"]["kraken"]["db"]
    output:
        "kraken/{unit}.tsv"
    conda:
        "../envs/kraken.yaml"
    log:
        "logs/kraken/{unit}.log"
    benchmark:
        "benchmarks/kraken/{unit}.txt"
    threads: 64
    params:
        gzip=lambda wildcards, input: "--gzip-compressed" if input.reads[0].endswith(".gz") else ""
    shell:
        """
        if [[ -s {input.reads[0]} ]]
        then
            kraken --fastq-input --paired {params.gzip} \
            --threads {threads} --db {input.db} \
            {input.reads} > {output} 2> {log}
        else
            touch {output}
        fi
        """


rule kraken_report:
    input:
        tsv="kraken/{unit}.tsv",
        db=config["params"]["kraken"]["db"]
    output:
        report("tables/{unit}.classification.tsv",
               caption="../report/kraken.rst", category="QC")
    conda:
        "../envs/kraken.yaml"
    log:
        "logs/kraken-report/{unit}.log"
    benchmark:
        "benchmarks/kraken/{unit}_report.txt"
    shell:
        "kraken-report --db {input.db} {input.tsv} > {output} 2> {log}"


rule calc_tree_colormap:
    input:
        "tables/{unit}.classification.tsv"
    output:
        "kraken/{unit}.colormap.pickle"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/calc-colormap.py"


rule extract_classification_tree:
    input:
        classification="tables/{unit}.classification.tsv",
        colormap="kraken/{unit}.colormap.pickle"
    output:
        "kraken/{unit}.classification.dot"
    benchmark:
        "benchmarks/kraken/{unit}_classification_tree.txt"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-classification.py"


rule plot_kmer_mapping:
    input:
        mapping="kraken/{unit}.tsv",
        colormap="kraken/{unit}.colormap.pickle"
    output:
        report("plots/{unit}.kmer-mapping.svg",
               caption="../report/kmer-mapping.rst", category="QC")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-kmer-mapping.py"


rule plot_classification_tree:
    input:
        "kraken/{unit}.classification.dot"
    output:
        report("plots/{unit}.classification.svg",
               caption="../report/classification-tree.rst",
               category="QC")
    conda:
        "../envs/eval.yaml"
    shell:
        "dot -Tsvg {input} "
        "-Grankdir=TB -Nshape=box -Nstyle=rounded -Nfontname=sans "
        "-Nfontsize=10 -Npenwidth=2 -Epenwidth=2 -Ecolor=grey -Nbgcolor=white "  # -Ncolor='{params.color}'"
        "> {output}"
