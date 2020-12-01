import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


rule download_minikraken_db:
    input:
        HTTP.remote(
            config["params"]["kraken"]["db_link"],
        )
    output:
        directory(config["params"]["kraken"]["db"]),
    params:
        minikraken_target=lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/python.yaml"
    log:
        "logs/kraken/db_download.log"
    shell:
        "tar -xvf {input} -C {params.minikraken_target} &> {log}"


rule kraken:
    input:
        reads=lambda w: units.loc[w.unit, ["fq1", "fq2"]],
        db=config["params"]["kraken"]["db"]
    output:
        "resources/kraken/{unit}.tsv"
    log:
        "logs/kraken/{unit}.log"
    conda:
        "../envs/kraken.yaml"
    benchmark:
        "benchmarks/kraken/{unit}.txt"
    threads:
        min(workflow.cores, 64)
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
        tsv="resources/kraken/{unit}.tsv",
        db=config["params"]["kraken"]["db"]
    output:
        report(
            "analysis/kraken/tables/{unit}.classification.tsv",
            caption="../report/kraken.rst",
            category="QC",
            subcategory="Kraken",
        )
    log:
        "logs/kraken-report/{unit}.log"
    conda:
        "../envs/kraken.yaml"
    benchmark:
        "benchmarks/kraken/{unit}_report.txt"
    shell:
        "kraken-report --db {input.db} {input.tsv} > {output} 2> {log}"


rule calc_tree_colormap:
    input:
        "analysis/kraken/tables/{unit}.classification.tsv"
    output:
        "analysis/kraken/{unit}.colormap.pickle"
    log:
        "logs/calc_tree_colormap/{unit}.log"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/calc-colormap.py"


rule extract_classification_tree:
    input:
        classification="analysis/kraken/tables/{unit}.classification.tsv",
        colormap="analysis/kraken/{unit}.colormap.pickle"
    output:
        "analysis/kraken/{unit}.classification.dot"
    benchmark:
        "benchmarks/kraken/{unit}_classification_tree.txt"
    log:
        "logs/extract_classification_tree/{unit}.log"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-classification.py"


rule plot_kmer_mapping:
    input:
        mapping="resources/kraken/{unit}.tsv",
        colormap="analysis/kraken/{unit}.colormap.pickle"
    output:
        report(
            "results/plots/{unit}.kmer-mapping.svg",
            caption="../report/kmer-mapping.rst",
            category="QC",
            subcategory="Kraken",
        )
    log:
        "logs/plot_kerm_mapping/{unit}.kmer-mapping.log"
    conda:
        "../envs/eval.yaml"
    script:
        "../scripts/plot-kmer-mapping.py"


rule plot_classification_tree:
    input:
        "analysis/kraken/{unit}.classification.dot"
    output:
        report(
            "results/plots/{unit}.classification.svg",
            caption="../report/classification-tree.rst",
            category="QC",
            subcategory="Kraken",
        )
    log:
        "logs/plot_classification_tree/{unit}.classification.log"
    conda:
        "../envs/eval.yaml"
    shell:
        "dot -Tsvg {input} "
        "-Grankdir=TB -Nshape=box -Nstyle=rounded -Nfontname=sans "
        "-Nfontsize=10 -Npenwidth=2 -Epenwidth=2 -Ecolor=grey -Nbgcolor=white "  # -Ncolor='{params.color}'"
        "> {output} 2> {log}"
