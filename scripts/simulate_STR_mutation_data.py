#!/usr/bin/env python3
import argparse

import numpy as np
import pandas as pd

def parse_cla():
   parser = argparse.ArgumentParser()

   parser.add_argument(
       "-r", "--str_regions", type=str, required=True, 
       help="File with STR coordinates and reference allele length to simulate reads for."
   )
   parser.add_argument(
       "-o", "--output", type=str, required=True, 
       help="Filepath to use for writing the simulated data."
   )
   parser.add_argument(
       "-n", "--n_patients", type=int, required=True, 
       help="The number of patients to simulate data for."
   )
   parser.add_argument(
       "-p", "--p_mut", type=float, required=True, 
       help="Probability for each STR locus to be mutated"
   )
   parser.add_argument(
       "-b", "--direction_bias", type=float, required=False, default=0.5,
       help="Proportion of STR mutations that result in a deletion (default=0.5, which means insertions and deletions are equally likely)."
   )
   parser.add_argument(
       "-z", "--stepsize", type=float, required=False, default=0.555,
       help="Mutation stepsizes are drawn from a geometric distribution,\
          --stepsize is used to set the probability of succes for the geometric \
            distribution that is used. (default=0.555, the observed value for STR mutations in MSS patients)"
   )
   parser.add_argument(
       "-s", "--random_seed", type=int, required=False, default=None,
       help="Seed to use for random number generation."
   )

   return parser.parse_args()

def main():
    args = parse_cla()

    df_loci = pd.read_csv(args.str_regions)

    patients = np.repeat(
        np.array([i for i in range(1, args.n_patients + 1)]), 
        len(df_loci)
    )
    df_str_lengths = pd.DataFrame(
        {
            "patient": patients,
            "tmp_id": np.tile(df_loci["tmp_id"], args.n_patients),
            "period": np.tile(df_loci["period"], args.n_patients),
            "ref": np.tile(df_loci["ref"], args.n_patients),
        }        
    )

    # for convenience' sake, all STR loci are homozygous WT in healthy samples
    df_str_lengths["allele_a_healthy"] = df_str_lengths["ref"]
    df_str_lengths["allele_b_healthy"] = df_str_lengths["ref"]
    df_str_lengths["allele_a_tumor"] = df_str_lengths["ref"]
    df_str_lengths["allele_b_tumor"] = df_str_lengths["ref"]    

    df_str_lengths = df_str_lengths.melt(
        id_vars=("patient", "tmp_id", "period", "ref"),
        value_vars=("allele_a_healthy", "allele_b_healthy", "allele_a_tumor", "allele_b_tumor"),
        var_name="allele",
        value_name="allele_length"
    )
    tumor_idx = df_str_lengths["allele"].str.endswith("_tumor")

    rng = np.random.default_rng(seed = args.random_seed)
    # 0.0142 for MSS, 0.0518 for MSI
    mutated = rng.choice([True, False], size = len(df_str_lengths[tumor_idx]), p=[args.p_mut, 1 - args.p_mut])
    allele_diff = np.zeros_like(mutated, dtype=int)
    direction = rng.choice([-1, 1], size = mutated.sum(), p=[args.direction_bias, 1 - args.direction_bias])

    # 0.555 MSS, 0.354 MSI
    allele_diff[mutated] = rng.geometric(args.stepsize, size = len(allele_diff[mutated])) * direction
    df_str_lengths.loc[tumor_idx, "allele_length"] += allele_diff

    # STR allele lengths should never be < 1
    df_str_lengths.loc[df_str_lengths["allele_length"] < 1, "allele_length"] = 1

    df_str_lengths = df_str_lengths.pivot(
        index=("patient", "tmp_id", "period", "ref" ),
        columns="allele",
        values="allele_length"
    ).reset_index().sort_values(["patient", "tmp_id"])

    df_str_lengths.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
