#!/usr/bin/env python3
import argparse

import pandas as pd

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mss", type=str, required=True,
        help="Filepath to dummy MSS data."
    )
    parser.add_argument(
        "--msi", type=str, required=True,
        help="Filepath to dummy MSI data"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, 
        help="Filepath to use for writing the merged data."
    )

    return parser.parse_args()

def main():
    args = cla_parser()

    df_mss = pd.read_csv(args.mss)
    df_mss["MSI"] = "MSS"
    n_mss = df_mss["patient"].nunique()
    
    df_msi = pd.read_csv(args.msi)
    df_msi["MSI"] = "MSI"
    df_msi.loc[:, "patient"] += n_mss

    df_merged = pd.concat([df_mss, df_msi])

    df_merged.to_csv(args.output, index=False)




if __name__ == "__main__":
    main()
