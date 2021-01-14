# Read in individuals and units
individuals = pd.read_csv(config["individuals"], sep="\t", dtype=str,).set_index(
    "id", drop=False
)
individuals["hash"] = np.arange(len(individuals))
units = pd.read_csv(config["units"], sep="\t", dtype=str,).set_index("unit", drop=False)


# Specify where to find the kraken database.
# Assemble the paths for the kraken targets, if kraken is used.
kraken_db = config["params"]["kraken"].get("db")
kraken_targets = []
if kraken_db:
    kraken_targets = expand(
        [
            "results/plots/{unit}.kmer-mapping.svg",
            "results/plots/{unit}.classification.svg",
        ],
        unit=units.unit,
    )


def param_str(s):
    """Assemble a parameter string in the format n=6.M=5.m=3
    from a stacks params dictionary in the config file.
    """
    return f"n={s['max_locus_mm']}.M={s['max_individual_mm']}.m={s['min_reads']}"


# Generate output files using the file extensions specified in the
# config file.
parameter_sets = [param_str(param_set) for param_set in config["params"]["stacks"]]


get_outdir = lambda w, output: os.path.dirname(output[0])


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
            print(
                f"Invalid output type {f_type} for populations. Should be one of {pop_file_suffixes.keys()}.",
                file=sys.stderr,
            )
    if suffixes:
        return suffixes
    else:
        print(f"No valid output files specified for populations.", file=sys.stderr)
        sys.exit(1)


def trim_input():
    """Depending on the desired mode, this can either be p5 reads only,
    p5 and p7 reads merged into one ('single end') read, or the p5 and p7
    files concatenated into one file.
    """
    mode = config["reads"]["mode"]
    if mode == "p5_only":
        return "analysis/extracted/{individual}.1.fq.gz"
    elif mode == "merged":
        return "analysis/merged/{individual}.fq.gz"
    elif mode == "concatenated":
        return "analysis/concatenated/{individual}.fq.gz"
    else:
        raise ValueError(
            f"Invalid mode: {mode}. Should be 'p5_only', 'merged', or 'concatenated'"
        )


def fmt_ustacks_input(wildcards, input):
    return ["-s {}".format(f[: -len(".tags.tsv.gz")]) for f in input.ustacks]


def get_cargopath(w, input):
    """Check if called from within the .test folder to get the right path.
    """
    if ".test" in os.path.abspath(input.tsv):
        return "../../workflow/scripts/stacks_analyzer/Cargo.toml"
    else:
        return "workflow/scripts/stacks_analyzer/Cargo.toml"


def get_all_parameter_sets():
    all_sets = []
    for parameters in config["params"]["stacks"]:
        all_sets.append(
            f"n={parameters['max_locus_mm']}.M={parameters['max_individual_mm']}.m={parameters['min_reads']}"
        )
    return all_sets
