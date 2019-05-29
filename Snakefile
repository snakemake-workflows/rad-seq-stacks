import pandas as pd
import numpy as np
import zlib
import sys

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

def pop_suffixes():
    """Map input file types of the stacks populations script
    to the file suffixes they generate.
    """
    pop_file_suffixes = {
        "fasta": ["samples-raw.fa"],
        "genepop": ["snps.genepop", "haps.genepop"],
        "vcf": ["snps.vcf", "haps.vcf"],
        "phylip": ["fixed.phylip", "fixed.phylip.log"],
    }
    suffixes = []
    for f_type in config["params"]["populations"]["output_types"]:
        try:
            suffixes.extend(pop_file_suffixes[f_type])
        except KeyError:
            print(f"Invalid output type {f_type} for populations. Should be one of {pop_file_suffixes.keys()}.", file=sys.stderr)
    if suffixes:
        return suffixes
    else:
        print(f"No valid output files specified for populations.", file=sys.stderr)
        sys.exit(1)


rule all:
    input:
        expand("calls/n={p[max_locus_mm]}.M={p[max_individual_mm]}.m={p[min_reads]}.populations.{e}",
               p=config["params"]["stacks"],
               e=pop_suffixes(),
        ),
        kraken_targets


include: "rules/common.smk"
include: "rules/preprocessing.smk"
include: "rules/stacks.smk"
include: "rules/kraken.smk"
