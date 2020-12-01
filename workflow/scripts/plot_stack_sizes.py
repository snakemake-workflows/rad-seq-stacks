""" Visualize stack size (reads per locus) distribution.
"""
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("agg")
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


def read_dat_file(dat_file_path, threshold):
    """Read dat files into a list. Split off everything over
    threshold into a separate list.
    """
    params = dat_file_path.split("/")[-2]
    name = os.path.splitext(os.path.basename(dat_file_path))[0]
    data = []
    outliers = []
    with open(dat_file_path, "r") as dat_file:
        for line in dat_file:
            count = int(line.strip())
            if count <= threshold:
                data.append(int(count))
            else:
                outliers.append(int(count))
    return data, params, name


def compare_parameters():
    """Plot distribution of Locus sizes for several parameter sets."""
    # thin lines
    mpl.rcParams['axes.linewidth'] = 0.5

    # parse input data
    df = None
    counts_df = None
    all_params = []
    cluster_counts = []
    threshold = snakemake.params.threshold

    for data_file_path in snakemake.input.dats:
        data, params, name = read_dat_file(data_file_path, threshold)
        file_df = pd.DataFrame(
            {
                "count": data,
                "params": [params for _ in enumerate(data)],
                "name": [name for _ in enumerate(data)],
            },
        )

        locus_count_df = pd.DataFrame(
            {
                "nr_loci": [len(data)],
                "params": [params],
                "name": [name],
            },
        )

        if df is None:
            df = file_df
        else:
            df = df.append(file_df, ignore_index=True)

        if counts_df is None:
            counts_df = locus_count_df
        else:
            counts_df = counts_df.append(locus_count_df, ignore_index=True)

        all_params.append(params)
        cluster_counts.append(len(data))

    df.to_csv(snakemake.output.sizes_dataframe)
    counts_df.to_csv(snakemake.output.counts_dataframe)

    scale = snakemake.params.scale
    fig, ax_violin = plt.subplots(
        figsize=(12*scale, 6*len(set(all_params))*scale)
    )

    thresholded_df = df[df["count"] < threshold]

    sns.violinplot(
        data=thresholded_df,
        x="count",
        y="params",
        ax=ax_violin,
        bw=snakemake.params.violin_plot_bw,
    )

    plt.tight_layout()
    plt.savefig(snakemake.output.violin_pdf, format="pdf", dpi=300)

    plt.clf()
    fig, ax_violin = plt.subplots(
        figsize=(8*scale, 8*scale)
    )
    sns.scatterplot(x="params", y="nr_loci", hue="name", data=counts_df)
    plt.savefig(snakemake.output.scatter_pdf, format="pdf", dpi=300)


if __name__ == "__main__":
    compare_parameters()
