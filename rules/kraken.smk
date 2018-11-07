rule kraken:
    input:
        fq1="trimmed-adapter/{individual}.1.fq.gz",
        fq2="trimmed-adapter/{individual}.2.fq.gz",
        db=config["params"]["kraken"]["db"]
    output:
        "kraken/{individual}.tsv"
    conda:
        "../envs/kraken.yaml"
    log:
        "logs/kraken/{individual}.log"
    threads: 64
    shell:
        """
        if [ -s {input.fq1} ]
        then
            kraken --threads {threads} --db {input.db} {input.fq1} {input.fq2} > {output} 2> {log}
        else
            touch {output}
        fi
        """


rule kraken_report:
    input:
        tsv="kraken/{individual}.tsv",
        db=config["params"]["kraken"]["db"]
    output:
        report("tables/{individual}.classification.tsv", caption="../report/kraken.rst", category="classification")
    conda:
        "../envs/kraken.yaml"
    log:
        "logs/kraken-report/{individual}.log"
    shell:
        "kraken-report --db {input.db} {input.tsv} > {output} 2> {log}"


rule calc_tree_colormap:
    input:
        "tables/{individual}.classification.tsv"
    output:
        "kraken/{individual}.colormap.pickle"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/calc-colormap.py"


rule extract_classification_tree:
    input:
        classification="tables/{individual}.classification.tsv",
        colormap="kraken/{individual}.colormap.pickle"
    output:
        "kraken/{individual}.classification.dot"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-classification.py"


rule plot_kmer_mapping:
    input:
        mapping="kraken/{individual}.tsv",
        colormap="kraken/{individual}.colormap.pickle"
    output:
        report("plots/{individual}.kmer-mapping.svg", caption="../report/kmer-mapping.rst", category="K-mer mapping")
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-kmer-mapping.py"


rule plot_classification_tree:
    input:
        "kraken/{individual}.classification.dot"
    output:
        report("plots/{individual}.classification.svg", caption="../report/classification-tree.rst", category="classification")
    conda:
        "../envs/eval.yaml"
    shell:
        "dot -Tsvg {input} "
        "-Grankdir=TB -Nshape=box -Nstyle=rounded -Nfontname=sans "
        "-Nfontsize=10 -Npenwidth=2 -Epenwidth=2 -Ecolor=grey -Nbgcolor=white " # -Ncolor='{params.color}'"
        "> {output}"
