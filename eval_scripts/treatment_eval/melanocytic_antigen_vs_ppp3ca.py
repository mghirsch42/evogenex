import pandas as pd
import numpy as np
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests

pearson = True
diff_exp = False

df = pd.read_csv("data/mouse_treatment/Melanocytic_antigen_in_R_vs_NR.csv", header=[0,1], index_col=0)
df = df.stack().stack().reset_index()
df.columns = ["gene", "response", "mouse", "expr"]
df["response2"] = ["responder" if df.loc[i, "response"] == "R" or df.loc[i, "response"] == "SD" 
                        else "non-responder"
                        for i in range(len(df))]

# Normalize
scaler_df = pd.read_csv("data/mouse_treatment/neutral_scaler_values.csv")
scaler_df["mouse"] = scaler_df["mouse"].astype(str)
scaled_expr = []
for idx, row in df.iterrows():
    scale_val = scaler_df[scaler_df["mouse"] == row["mouse"]]["scale_val"].values[0]
    scaled_expr.append(row["expr"] / scale_val)
df["expr"] = scaled_expr

# Pearson correlations with Ppp3ca
if pearson:
    ppp3ca = df[df["gene"] == "Ppp3ca"]
    df = df[df["gene"] != "Ppp3ca"]
    genes = []
    pearson_stats = []
    pearson_pvalues = []
    for gene in df["gene"].unique():
        pearson_res = stats.pearsonr(ppp3ca["expr"].values, df[df["gene"] == gene]["expr"].values)
        genes.append(gene)
        pearson_stats.append(round(pearson_res.statistic, 5))
        pearson_pvalues.append(round(pearson_res.pvalue, 5))

    pearson_results = pd.DataFrame()
    pearson_results["gene"] = genes
    pearson_results["pearson_stat"] = pearson_stats
    pearson_results["pearson_pvalue"] = pearson_pvalues
    print(pearson_results)

    # pearson_results.to_csv("results/mouse_treatment/wnt/melanocytic_antigen_vs_ppp3ca.csv", index=False)

# Load the t-test resutls for these individual genes
if diff_exp:
    diff_exp_df = pd.read_csv("results/mouse_treatment/norm_ttest.csv")
    genes = df["gene"].unique()
    t_df = diff_exp_df[diff_exp_df["gene_name"].isin(genes)]
    print(t_df)
    t_df.to_csv("results/mouse_treatment/wnt/melanocytic_antigen_ttest_data.csv", index=False)
