import pandas as pd
import pysam


def main(args):
    bam = pysam.AlignmentFile(args.bam_file)

    bl_df = pd.read_csv(args.black_list_file, header=None, sep="\t")

    bl_df.columns = "chrom", "beg", "end", "comment"

    chroms = args.chromosomes

    if chroms is None:
        chroms = bam.references

    df = []

    for chrom in chroms:
        df.extend(
            get_chromosome_coverage(
                bam, bl_df[bl_df["chrom"] == chrom], chrom, bin_size=args.bin_size, min_mqual=args.min_mqual
            )
        )

    df = pd.DataFrame(df)

    df.to_csv(args.out_file, compression="gzip", index=False, sep="\t")


def get_chromosome_coverage(bam, bl_df, chrom, bin_size=int(5e5), min_mqual=30):
    df = []

    chrom_len = bam.get_reference_length(chrom)

    for i in range(0, chrom_len, bin_size):
        beg = i

        end = i + bin_size

        valid = True

        for _, row in bl_df.iterrows():
            if (row["beg"] <= beg <= row["end"]) or (row["beg"] <= end <= row["end"]):
                valid = False

        if valid:
            count = bam.count(
                chrom,
                beg,
                end,
                read_callback=lambda x: (check_read(x, min_mqual=min_mqual)),
            )

            df.append({"chrom": chrom, "start": beg, "end": end, "reads": count})

    return df


def check_read(read, min_mqual=30):
    valid = True

    if read.mapping_quality < min_mqual:
        valid = False

    elif read.is_duplicate:
        valid = False

    elif read.is_unmapped:
        valid = False

    elif read.is_qcfail:
        valid = False

    elif read.is_secondary:
        valid = False

    return valid


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bam-file", required=True)

    parser.add_argument("-k", "--black-list-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-c", "--chromosomes", nargs="+", type=str, default=None)

    parser.add_argument("-q", "--min-mqual", type=int, default=30)

    parser.add_argument("-s", "--bin-size", type=int, default=int(5e5))

    cli_args = parser.parse_args()

    main(cli_args)
