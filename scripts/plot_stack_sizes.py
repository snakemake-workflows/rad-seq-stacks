""" Visualize stack size (reads per locus) distribution.

Usage:
    python plot_stack_sizes.py stacks_sizes file_1.dat file_2.dat file_3.dat
"""
import click
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd


def tuftefy(ax):
    """Remove spines and tick position markers to reduce ink."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color('grey')
    ax.grid(color="w", alpha=0.7)
    ax.get_yaxis().grid(True)
    ax.get_xaxis().grid(False)


def dat_file_to_array(dat_file, threshold):
    """Read dat files into a list. Split off everything over
    threshold into a separate list.
    """
    data = []
    outliers = []
    for line in dat_file:
        count = int(line.strip())
        if count <= threshold:
            data.append(int(count))
        else:
            outliers.append(int(count))
    return data, outliers


@click.group()
def cli():
    """Plot histogram and violin plot of stack size distributions."""
    # Note: This has to stand above the other click code.
    pass


@cli.command("stacks_sizes", short_help="Visualize distribution of stack "
             "sizes extracted from stacks tsv files.")
@click.argument("data_file", nargs=1, type=click.File(mode='r'))
@click.option("--output-path", type=click.Path(exists=False, dir_okay=False))
@click.option("--threshold", type=int, default=500, help="Threshold to cut off stacks sizes in the histogram")
@click.option("--ylim", type=int, default=None, help="Force the plot to have a fixed upper value for the y-axis.")
def stack_size_distribution(data_file, output_path, threshold, ylim):
    """Distribution of one or more
    """
    # thin lines
    mpl.rcParams['axes.linewidth'] = 0.5
    # parse input data
    data, outliers = dat_file_to_array(data_file, threshold)
    # Set up plots to contain the histogram (left) and violin plots (right)
    # in approx 4:3 ratio
    fig, (ax_hist, ax_violin) = plt.subplots(nrows=1, ncols=2, figsize=(18, 8),
                                             gridspec_kw={
                                                 'width_ratios': [3, 1]
                                             })
    ax_hist.set_title(os.path.basename(data_file.name))

    # plot main histogram
    ax_hist.hist(data, bins=min(threshold, 1000), alpha=0.7, align="mid",
                 label=f"{data_file.name}")
    # denote threshold line
    ax_hist.axvline(threshold, color="gray",
                    linewidth=1).set_linestyle("dashed")
    tuftefy(ax_hist)

    # for the same y limit on all plots, if the parameter has been provided
    if ylim is not None:
        ax_hist.set_ylim((0, ylim))
    ax_hist.set_xlabel(f"Stack sizes (truncated at {threshold})")
    ax_hist.set_ylabel(f"Number of observed stacks")
    ax_hist.legend()

    d = pd.DataFrame(
        {
            "count": data + outliers,
            "params": ["" for _ in data + outliers],  # used to show at one spot
        }
    )
    # logarithmic violin plot as secondary view on the data
    sns.violinplot(data=d, x="params", y="count", ax=ax_violin, inner="quartile", bw=0.2)
    ax_violin.set_yscale('log')

    # force upper y-axis limit to make plots comparable
    ax_violin.set_ylim((1, 10**5))
    ax_violin.set_ylabel(f"Stack sizes (log; not truncated)")

    if output_path:
        print(f"Writing plot to: {output_path}")
    else:
        output_path = os.path.splitext(data_file.name)[0] + ".pdf"
    plt.tight_layout()
    plt.savefig(output_path, format="pdf", dpi=70)


if __name__ == "__main__":
    cli()
