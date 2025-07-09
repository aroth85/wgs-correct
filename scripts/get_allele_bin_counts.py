import pandas as pd
import pysam


def main(args):
    bam = pysam.AlignmentFile(args.bam_file, "r")

    bl_df = pd.read_csv(args.black_list_file, header=None, sep="\t")

    bl_df.columns = "chrom", "beg", "end", "comment"

    df = pd.read_csv(args.in_file, converters={"chrom": str, "coord": int}, sep="\t")

    out_df = []

    # TODO: Replace O(n^2) with proper O(n) calc
    for chrom in args.chromosomes:
        print(chrom)

        bl_chrom_df = bl_df[bl_df["chrom"] == chrom]

        chrom_df = df[df["chrom"] == chrom]

        chrom_len = bam.get_reference_length(chrom)

        for beg in range(0, chrom_len, args.bin_size):
            end = beg + args.bin_size

            valid = True

            for _, row in bl_chrom_df.iterrows():
                if (row["beg"] <= beg <= row["end"]) or (row["beg"] <= end <= row["end"]):
                    valid = False
                    break

            if valid:
                bin_df = chrom_df[chrom_df["coord"].between(beg, end, inclusive="both")]

                row = {
                    "chrom": chrom,
                    "start": beg,
                    "end": end,
                    "ref_count": bin_df["ref_count"].sum(),
                    "alt_count": bin_df["alt_count"].sum(),
                    "allele_0_count": bin_df["allele_0_count"].sum(),
                    "allele_1_count": bin_df["allele_1_count"].sum(),
                }

                out_df.append(row)

    out_df = pd.DataFrame(out_df)

    if args.sample_id is not None:
        out_df.insert(0, "sample", args.sample_id)

    out_df.to_csv(args.out_file, compression="gzip", index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bam-file", required=True)

    parser.add_argument("-c", "--chromosomes", nargs="+", required=True)

    parser.add_argument("-i", "--in-file", required=True)

    parser.add_argument("-k", "--black-list-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-s", "--bin-size", required=True, type=int)

    parser.add_argument("--sample-id", default=None)

    cli_args = parser.parse_args()

    main(cli_args)
