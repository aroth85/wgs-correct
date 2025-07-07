import matplotlib.pyplot as pp
import numpy as np
import pandas as pd
import seaborn as sb


def main(args):
    df = pd.read_csv(args.in_file, converters={"chrom": str}, sep="\t")

    df["baf"] = df["alt_count"] / (df["ref_count"] + df["alt_count"])

    df["baf_phased"] = df["allele_1_count"] / (df["allele_0_count"] + df["allele_1_count"])

    chroms = sort_chroms(df["chrom"].unique())

    chroms_size = df["chrom"].value_counts()

    width_ratios = [chroms_size[x] for x in chroms]

    plot_vals = ["baf", "baf_phased"]

    fig = pp.figure(figsize=(16, 4))

    grid = fig.add_gridspec(2, 1, hspace=0.1)

    for i, v in enumerate(plot_vals):
        sub_grid = grid[i].subgridspec(
            1, len(chroms), width_ratios=width_ratios, wspace=0.05
        )

        plot_by_chrom(df, chroms, fig, sub_grid, title=v, y_col=v)

        grid.tight_layout(fig)

    fig.savefig(args.out_file, dpi=120, bbox_inches="tight")


def plot_by_chrom(df, chroms, fig, grid, title=None, y_col="baf"):
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

        m = df[y_col].mean()

        s = df[y_col].std()

        ax.set_ylim(max(m - 3 * s, 0), min(m + 3 * s, 1))

    if title is not None:
        ax = fig.add_subplot(grid[:])

        ax.axis("off")

        ax.set_title(title)


def sort_chroms(chroms):
    numeric = []

    string = []

    if chroms[0].startswith("chr"):
        chr_prefix = True

    else:
        chr_prefix = False

    for c in chroms:
        if chr_prefix:
            c = c.replace("chr", "")

        try:
            numeric.append(int(c))

        except ValueError:
            string.append(c)

    chroms = [str(x) for x in sorted(numeric)] + list(sorted(string))

    if chr_prefix:
        chroms = ["chr{}".format(x) for x in chroms]

    return chroms


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
