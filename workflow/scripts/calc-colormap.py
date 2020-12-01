import pickle
import common
import sys

with open(snakemake.log[0], "w") as log_file:
    # redirect all output top the logfile as suggested
    # here: https://stackoverflow.com/a/64114743/2862719
    # sys.stderr = sys.stdout = log_file
    full_tree = common.load_classification_tree(
        snakemake.input[0],
        min_percentage=0.0,
    )
    cmap = common.TreeCMap(full_tree)
    with open(snakemake.output[0], "wb") as out:
        pickle.dump(cmap, out)
