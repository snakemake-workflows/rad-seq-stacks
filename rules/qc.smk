def get_cargopath(w, input):
    """Check if called from within the .test folder to get the right path.
    """
    if ".test" in os.path.abspath(input.tsv):
        return "../scripts/stacks_analyzer/Cargo.toml"
    else:
        return "scripts/stacks_analyzer/Cargo.toml"

def get_plot_script(w, input):
    prefix = "../" if ".test" in os.path.abspath(input.dats) else ""
    return prefix + "scripts/plot_stack_sizes.py"


# compute the lengths of stacks loci by counting lines in the tags.tsv files
rule count:
    input:
        tsv="stacks/{parameter_set}/{individual}.tags.tsv.gz"
    output:
        dat="counts/{parameter_set}/{individual}.dat"
    params:
        cargo_path=get_cargopath
    conda:
        "../envs/rust.yaml"
    shell:
        """
        cargo run --release --manifest-path={params.cargo_path} -- count-locus-sizes {input.tsv} --dat {output.dat}
        """


# plot a histogram and a violin plot of stacks sizes for one individual
rule plot:
    input:
        dats="counts/{parameter_set}/{individual}.dat",
    output:
        pdf="plots/{parameter_set}/{individual}.pdf"
    conda:
        "../envs/plot_stacks_dist.yaml"
    params:
        plot_script=get_plot_script
    shell:
        """
        python {params.plot_script} stacks_sizes {input.dats} --output-path {output.pdf}
        """


# concatenate all plots for one parameter set into one pdf document
rule assemble_report:
    input:
        plots=expand("plots/{{parameter_set}}/{individual}.pdf",
                     individual=individuals.id,
    )
    output:
        pdf="plots/stacks_size_distribution_{parameter_set}.pdf"
    conda:
        "../envs/pdfunite.yaml"
    shell:
        "pdfunite {input.plots} {output.pdf}"
