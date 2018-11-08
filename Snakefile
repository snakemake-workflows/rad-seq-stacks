import pandas as pd
import numpy as np
import zlib

configfile: "config.yaml"


individuals = pd.read_table("individuals.tsv", dtype=str).set_index("id", drop=False)
individuals["hash"] = np.arange(len(individuals)) # individuals.id.str.encode("utf-8").apply(zlib.crc32)
units = pd.read_table("units.tsv", dtype=str).set_index("id", drop=False)


rule all:
    input:
        expand("calls/n={p[max_locus_mm]}.M={p[max_individual_mm]}.m={p[min_reads]}/populations.snps.vcf",
               p=config["params"]["stacks"]),
        expand("validation/n={p[max_locus_mm]}.M={p[max_individual_mm]}.m={p[min_reads]}/validation.txt",
               p=config["params"]["stacks"]),
        expand(["plots/{unit}.classification.svg",
                "plots/{unit}.kmer-mapping.svg"],
               unit=units.id)


include: "rules/common.smk"
include: "rules/preprocessing.smk"
include: "rules/stacks.smk"
include: "rules/kraken.smk"
include: "rules/validation.smk"
