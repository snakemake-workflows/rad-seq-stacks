import pandas as pd

configfile: "config.yaml"


individuals = pd.read_table("individuals.tsv").set_index("id", drop=False)
units = pd.read_table("units.tsv").set_index("id", drop=False)


rule all:
    input:
        expand("gstacks/n={p[max_locus_mm]}.M={p[max_individual_mm]}.m={p[min_reads]}/calls.vcf",
               p=config["params"]["stacks"])


include: "rules/preprocessing.smk"
include: "rules/stacks.smk"
