#!/usr/bin/env python3
import argparse

import pandas as pd

def parse_cla():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-r", "--str_regions", type=str, required=True, 
        help="File with STR coordinates and reference allele length to filter."
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, 
        help="Filepath to use for writing the filtered data."
    )
    parser.add_argument(
        "-n", "--n_regions", type=int, required=False, 
        help="Number of regions to keep after filtering."
    )
    parser.add_argument(
        "-g", "--n_genes", type=int, required=False, 
        help="Number of distinct genes to simulate"
    )

    return parser.parse_args()

def char_generator():
    char = ord('A')
    while True:
        yield chr(char)
        char += 1
        if char > 90:
            raise ValueError("character generator out of bounds!")

def main():
    args = parse_cla()    

    df_loci = pd.read_csv(args.str_regions, sep="\t")
    n_regions = args.n_regions
    if not n_regions:
        n_regions = len(df_loci)

    df_loci_filt = (
        df_loci
            .query("neighbour_type == 'no_neighbour' and not in_segdup")
            .filter(['tmp_id','chr', 'start', 'end', 'period', 'ref'], axis=1)
            .head(n_regions)
    )

    if args.n_genes:
        bins = pd.qcut(df_loci_filt.start, args.n_genes)
        gene_gen = char_generator()
        bin_gene_map = dict()
        for i in bins.drop_duplicates():
            bin_gene_map[i] = next(gene_gen)
        df_loci_filt = df_loci_filt.assign(gene = lambda x: [bin_gene_map[i] for i in bins])

    df_loci_filt.to_csv(args.output, index=False)



if __name__ == "__main__":
    main()
