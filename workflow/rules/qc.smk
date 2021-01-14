# compute the lengths of stacks loci by counting lines in the tags.tsv files
rule count:
    input:
        tsv="analysis/stacks/{parameter_set}/{individual}.tags.tsv.gz",
    output:
        dat="analysis/counts/{parameter_set}/{individual}.dat",
    params:
        cargo_path=get_cargopath,
    log:
        "logs/count/{parameter_set}/{individual}.log",
    conda:
        "../envs/rust.yaml"
    shell:
        """
        cargo run --release --manifest-path={params.cargo_path} -- count-locus-sizes {input.tsv} --dat {output.dat} &> {log}
        """


rule plot_comparison:
    input:
        dats=expand(
            "analysis/counts/{parameter_set}/{individual}.dat",
            parameter_set=get_all_parameter_sets(),
            individual=individuals.id,
        ),
    output:
        violin_pdf=report(
            "results/plots/distribution_comparison/stacks_size_distribution.pdf",
            caption="../report/size_distribution.rst",
            category="QC",
            subcategory="Workflow",
        ),
        scatter_pdf=report(
            "results/plots/distribution_comparison/stacks_counts.pdf",
            caption="../report/size_distribution.rst",
            category="QC",
            subcategory="Workflow",
        ),
        sizes_dataframe="results/plots/distribution_comparison/sizes_dataframe.csv",
        counts_dataframe="results/plots/distribution_comparison/counts_dataframe.csv",
    log:
        "logs/plot_comparison.log",
    conda:
        "../envs/plot_stacks_dist.yaml"
    params:
        threshold=config["qc"]["violin_plots"]["threshold"],
        violin_plot_bw=config["qc"]["violin_plots"]["kernel_bandwith"],
        scale=config["qc"]["violin_plots"]["scale"],
    script:
        "../scripts/plot_stack_sizes.py"


rule post_run_qc:
    input:
        trimmed_spacer=expand(
            "analysis/trimmed-spacer/{unit}.2.fq.gz", unit=units.unit,
        ),
        consensus_fq1=expand(
            "analysis/dedup/{unit}.consensus.1.fq.gz", unit=units.unit,
        ),
        consensus_fq2=expand(
            "analysis/dedup/{unit}.consensus.2.fq.gz", unit=units.unit,
        ),
        extracted_fq1=expand(
            "analysis/extracted/{individual}.1.fq.gz", individual=individuals.id
        ),
        extracted_fq2=expand(
            "analysis/extracted/{individual}.2.fq.gz", individual=individuals.id
        ),
        extract_logs=expand("logs/extract/{unit}.log", unit=units.unit,),
    output:
        qc_log=report(
            "results/qc.log",
            caption="../report/workflow_qc.rst",
            category="QC",
            subcategory="Workflow",
        ),
    log:
        "logs/post_run_qc/post_run_qc.log",
    params:
        dbr_suffix=config["umi"]["fixed_suffix"],
        p5_residue=config["restriction-enzyme"]["p5"]["residue-seq"],
        p7_residue=config["restriction-enzyme"]["p7"]["residue-seq"],
        fail_on_error=config["qc"]["fail_on_error"],
        samples_per_file=config["qc"]["samples_per_file"],
        warning_threshold=config["qc"]["extract"]["warning_threshold"],
        error_threshold=config["qc"]["extract"]["error_threshold"],
    conda:
        "../envs/plot_stacks_dist.yaml"
    script:
        "../scripts/quality_control.py"
