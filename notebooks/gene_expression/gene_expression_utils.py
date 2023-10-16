#!/usr/bin/env python3
import numpy as np
import pandas as pd
import statsmodels.api as sm

def gt_genexp_correlation(args):
    gene, locus, df_genotypes, df_genexp = args

    df_results = {
        "gene": [],
        "tmp_id": [],
        "coefficient": [],
        "pvalue_coef": [],
        "intercept": []
    }

    target_genotypes = df_genotypes.query(f"tmp_id == '{locus}'")[['patient', 'mean_gt']]
    target_genexp = df_genexp.query(f"gene == '{gene}'")[['patient', 'expression']]

    df_gt_expression = target_genexp.merge(target_genotypes, on="patient", how="inner")
    if df_gt_expression['mean_gt'].nunique() < 3:
        return pd.DataFrame(df_results)
    
    X = sm.add_constant(df_gt_expression['mean_gt'])
    Y = df_gt_expression['expression']    
    
    model = sm.OLS(Y, X)
    result = model.fit()
    
    df_results['gene'].append(gene)
    df_results['tmp_id'].append(locus)
    df_results['coefficient'].append(result.params.iloc[1])
    df_results['pvalue_coef'].append(result.pvalues.iloc[1])
    df_results['intercept'].append(result.params.iloc[0])

    return pd.DataFrame(df_results)

def gt_genexp_correlation_permutated(args):
    gene, locus, df_genotypes, df_genexp = args

    df_results = {
        "gene": [],
        "tmp_id": [],
        "coefficient": [],
        "pvalue_coef": [],
        "intercept": []
    }
    
    rng = np.random.default_rng(42)
    target_genotypes = df_genotypes.query(f"tmp_id == '{locus}'")[['patient', 'mean_gt']]
    target_genotypes = target_genotypes.assign(mean_gt = lambda x: rng.permutation(x.mean_gt.values))
    target_genexp = df_genexp.query(f"gene == '{gene}'")[['patient', 'expression']]

    df_gt_expression = target_genexp.merge(target_genotypes, on="patient", how="inner")
    if df_gt_expression['mean_gt'].nunique() < 3:
        return pd.DataFrame(df_results)
    
    X = sm.add_constant(df_gt_expression['mean_gt'])
    Y = df_gt_expression['expression']    
    
    model = sm.OLS(Y, X)
    result = model.fit()
    
    df_results['gene'].append(gene)
    df_results['tmp_id'].append(locus)
    df_results['coefficient'].append(result.params[1])
    df_results['pvalue_coef'].append(result.pvalues[1])
    df_results['intercept'].append(result.params[0])

    return pd.DataFrame(df_results)