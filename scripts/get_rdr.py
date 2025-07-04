import pandas as pd


def main(args):
    cell_df = pd.read_csv(args.in_file, converters={"chrom": str}, sep="\t")

    norm_df = pd.read_csv(args.normal_file, converters={"chrom": str}, sep="\t")

    norm_df = norm_df.rename(columns={"reads": "normal"})

    df = pd.merge(cell_df, norm_df, on=["chrom", "start", "end"])

    scale = df["normal"].sum() / df["reads"].sum()

    ratio = df["reads"] / df["normal"]

    df["rdr"] = ratio * scale

    df = df[["chrom", "start", "end", "normal", "reads", "rdr"]]

    df.to_csv(args.out_file, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-n", "--normal-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    cli_args = parser.parse_args()

    main(cli_args)
