#!/usr/bin/env python3
import argparse

import numpy as np
import pandas as pd

# Dictionary to specify which STR loci will be eSTRs. eSTRs are specified 
# as <STR locus>: <eSTR effect> 
# The eSTR effect is interpreted as 'fraction of expression change per unit
# difference from the reference allele'.
# 
# E.g., for an eSTR with ref length 12 and eSTR effect of 0.2, if the eSTR 
# length increases to 13, expression will be 20% higher
# if the eSTR contracts to length 10, expression will be 40% lower.
# ESTRS = {
#     'chr1_9691434': -0.2,
#     'chr1_1295260': 0.2,
#     'chr1_10117623': -0.2,
#     'chr1_6834638': 0.2,
#     'chr1_10662707': -0.2,
#     'chr1_1047113': 0.2,
# }
GENE_ESTR_MAP = {
    'A':	('chr1_1047113', 4, -0.2),
    'B':	('chr1_1295260', 5, 0.2),
    'J':	('chr1_6834638', 18, -0.2),
    'M':	('chr1_9691434', 5, 0.2),
    'N':	('chr1_10117623', 12, -0.2),
    'P':	('chr1_10662707', 5, 0.2),
}

def parse_cla():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-r", "--str_regions", type=str, required=True, 
        help="File with STR coordinates and reference allele length to filter."
    )
    parser.add_argument(
        "-g", "--genotypes", type=str, required=True,
        help="Filepath to dummy STR genotyping data."
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, 
        help="Filepath to use for writing the filtered data."
    )
    parser.add_argument(
        "-s", "--random_seed", type=int, required=False, default=None,
        help="Seed to use for random number generation."
    )

    return parser.parse_args()


def main():
    args = parse_cla()

    # Load STR regions and genotyping data
    df_strs = pd.read_csv(args.str_regions)        
    df_somatic_mutations = pd.read_csv(args.genotypes
        ).assign(mean_gt = lambda x: (x.allele_a_tumor + x.allele_b_tumor) / 2
            ).merge(df_strs[["tmp_id", "gene"]], on="tmp_id")
    
    n_patients = df_somatic_mutations["patient"].nunique()
    genes = df_strs["gene"].drop_duplicates().sort_values().to_numpy()

    rng = np.random.default_rng(seed = args.random_seed)
    mean_expressions = rng.integers(1, 500, size=len(genes))

    expression_matrix = np.array([mean_expressions] * n_patients, dtype=np.float64).reshape((n_patients, len(genes)))
    for i, gene in enumerate(genes):
        noise = ((rng.uniform(size=n_patients) * 2) - 1) * (0.05 * mean_expressions[i])
        expression_matrix[:, i] += noise
        try:
            eSTR = GENE_ESTR_MAP[gene]
            genotypes = df_somatic_mutations.loc[df_somatic_mutations["tmp_id"] == eSTR[0]].sort_values("patient")["mean_gt"].values
            genotypes -= eSTR[1]
            expression_matrix[:, i] += (expression_matrix[:, i] * (eSTR[2] * genotypes))
        except KeyError:
            continue

    df_gene_expression = (
        pd.DataFrame(expression_matrix, columns=genes)
            .melt(var_name='gene', value_name='expression', ignore_index=False)
            .reset_index(drop=False)
            .rename(columns = {'index': 'patient'})
    )
    df_gene_expression.loc[:, "patient"] += 1

    df_gene_expression.to_csv(args.output, index=False)




    ### OLD CODE STARTS HERE ###
    # # Get gene names, simulate mean expression for each (one value per gene)
    # genes = df_strs["gene"].drop_duplicates().sort_values().to_numpy()
    # rng = np.random.default_rng(seed = args.random_seed)
    # mean_expressions = rng.integers(1, 500, size=len(genes))

    # # Add gene expression to genotyping data, initially: mean expression value for each gene
    # df_mean_expressions = pd.DataFrame({"gene": genes, "expression": mean_expressions})
    # df_somatic_mutations = df_somatic_mutations.merge(df_mean_expressions, on="gene")
    # print(df_somatic_mutations.shape)

    # # Simulate some noise to add to expression data
    # noise = ((rng.uniform(size=len(df_somatic_mutations)) * 2) - 1) * (0.01 * df_somatic_mutations["expression"])
    # df_somatic_mutations.loc[:, "expression"] += noise

    # # Add eSTR effects to those loci designated as eSTRs by dictionary
    # df_somatic_mutations = (
    #     df_somatic_mutations
    #         .assign(
    #             coef = lambda x: [ESTRS[i] if i in ESTRS.keys() else 0 for i in x.tmp_id],
    #             diff_to_ref = lambda x: x.mean_gt - x.ref,
    #             expression = lambda x: x.expression + (x.expression * (x.coef * x.diff_to_ref)))
    #         .drop(["coef", "diff_to_ref"], axis=1)
    # )
    # print(df_somatic_mutations.shape)

    # df_gene_expression = df_somatic_mutations[["patient", "gene", "expression"]].groupby(["patient", "gene"], as_index=False).agg("mean")
    # df_gene_expression.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
