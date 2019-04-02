import pandas as pd
import numpy as np
import zlib

configfile: "config.yaml"


individuals = pd.read_csv("individuals.tsv", sep="\t", dtype=str).set_index("id", drop=False)
individuals["hash"] = np.arange(len(individuals)) # individuals.id.str.encode("utf-8").apply(zlib.crc32)
units = pd.read_csv("units.tsv", sep="\t", dtype=str).set_index("id", drop=False)

kraken_db = config["params"]["kraken"].get("db")
kraken_targets = []
if kraken_db:
    kraken_targets = expand(["plots/{unit}.kmer-mapping.svg",
                             "plots/{unit}.classification.svg"],
                            unit=units.id)


rule all:
    input:
        expand("calls/n={p[max_locus_mm]}.M={p[max_individual_mm]}.m={p[min_reads]}.populations.snps.vcf",
               p=config["params"]["stacks"]),
        kraken_targets


include: "rules/common.smk"
include: "rules/preprocessing.smk"
include: "rules/stacks.smk"
include: "rules/kraken.smk"
