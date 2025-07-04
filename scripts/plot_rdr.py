import matplotlib.pyplot as pp
import numpy as np
import pandas as pd
import seaborn as sb

from utils import sort_chroms


def main(args):
    df = pd.read_csv(args.in_file, converters={"chrom": str}, sep="\t")

    df["reads"] = np.log2(df["reads"] + 1) - np.log2((df["reads"] + 1).mean())

    chroms = sort_chroms(df["chrom"].unique())

    chroms_size = df["chrom"].value_counts()

    width_ratios = [chroms_size[x] for x in chroms]

    plot_vals = ["reads", "rdr", "rdr_cor"]

    fig = pp.figure(figsize=(16, 6))

    grid = fig.add_gridspec(3, 1, hspace=0.1)

    for i, v in enumerate(plot_vals):
        sub_grid = grid[i].subgridspec(
            1, len(chroms), width_ratios=width_ratios, wspace=0.05
        )

    plot_by_chrom(df, chroms, fig, sub_grid, title=v, y_col=v)

    grid.tight_layout(fig)

    fig.savefig(args.out_file, dpi=120, bbox_inches="tight")


def plot_by_chrom(df, chroms, fig, grid, title=None, y_col="rdr"):
    for i, chrom in enumerate(chroms):
        chrom_df = df[df["chrom"] == chrom]

        chrom_df = chrom_df.sort_values(by=["start"])

        num_bins = chrom_df.shape[0]

        chrom_df["idx"] = np.arange(num_bins)

        ax = fig.add_subplot(grid[0, i])

        ax.scatter(
            np.arange(num_bins),
            chrom_df[y_col],
            s=5,
        )

        sb.despine(ax=ax, offset=10)

        ax.spines["top"].set_visible(False)

        ax.spines["right"].set_visible(False)

        if i != 0:
            ax.spines["left"].set_visible(False)

            ax.set_yticks([])

            ax.set_yticklabels([])

        else:
            ax.tick_params(axis="x", which="major", labelsize=12)

        ax.set_xticks([num_bins / 2])

        ax.set_xticklabels([chrom.replace("chr", "")], fontsize=12)

        try:
            ax.set_ylim(df[y_col].quantile(0.01), df[y_col].quantile(0.99))

        except ValueError:
            continue

    if title is not None:
        ax = fig.add_subplot(grid[:])

        ax.axis("off")

        ax.set_title(title)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
