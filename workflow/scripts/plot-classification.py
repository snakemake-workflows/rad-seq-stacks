import pickle
from networkx.drawing.nx_agraph import write_dot
import common
import sys


with open(snakemake.log[0], "w") as log_file:
    # redirect all output top the logfile as suggested
    # here: https://stackoverflow.com/a/64114743/2862719
    sys.stderr = sys.stdout = log_file
    cmap = pickle.load(open(snakemake.input.colormap, "rb"))

    tree = common.load_classification_tree(
        snakemake.input.classification,
        min_percentage=1.0,
    )
    cmap.color_tree(tree)

    write_dot(tree, snakemake.output[0])
