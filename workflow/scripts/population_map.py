import sys

with open(snakemake.log[0], "w") as log_file:
    # redirect all output top the logfile as suggested
    # here: https://stackoverflow.com/a/64114743/2862719
    sys.stderr = sys.stdout = log_file

    individuals = snakemake.params.individuals
    d = individuals[["id"]]
    d["pop"] = 1
    d.to_csv(snakemake.output[0], index=False, header=None, sep="\t")
