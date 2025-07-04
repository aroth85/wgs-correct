import pysam
import pandas as pd


def main(args):
    bam_reader = pysam.AlignmentFile(args.bam_file, "rb")

    snp_reader = pysam.VariantFile(args.snps_file, "r")

    nuc_map = {"A": 0, "C": 1, "G": 2, "T": 3}

    df = []

    for chrom in args.chromosomes:
        for record in snp_reader.fetch(contig=chrom):
            # Skip indels
            if len(record.alleles[0]) != 1 or len(record.alleles[1]) != 1:
                continue

            sample_id = list(record.samples.keys())[0]

            sample_info = record.samples[sample_id]

            # Skip homozygous SNPs
            if sample_info["GT"] not in set([(0, 1), (1, 0)]):
                continue

            counts = bam_reader.count_coverage(
                record.chrom,
                start=record.pos - 1,
                stop=record.pos,
                quality_threshold=args.min_base_qual
            )

            ref = record.alleles[0]

            alt = record.alleles[1]

            allele_0 = record.alleles[sample_info["GT"][0]]

            allele_1 = record.alleles[sample_info["GT"][1]]

            if allele_0 not in nuc_map or allele_1 not in nuc_map:
                continue

            row = {
                "chrom": record.chrom,
                "coord": record.pos,
                "ref": record.alleles[0],
                "alt": record.alleles[1],
                "ref_count": counts[nuc_map[ref]][0],
                "alt_count": counts[nuc_map[alt]][0],
                "allele_0": allele_0,
                "allele_1": allele_1,
                "allele_0_count": counts[nuc_map[allele_0]][0],
                "allele_1_count": counts[nuc_map[allele_1]][0]
            }

            df.append(row)

    df = pd.DataFrame(df)

    df.to_csv(args.out_file, index=False, sep="\t")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--chromosomes", nargs="+", type=str, required=True)

    parser.add_argument("-b", "--bam-file", required=True)

    parser.add_argument("-s", "--snps-file", required=True)

    parser.add_argument("-o", "--out-file", required=True)

    parser.add_argument("-q", "--min-base-qual", type=int, default=30)

    cli_args = parser.parse_args()

    main(cli_args)
