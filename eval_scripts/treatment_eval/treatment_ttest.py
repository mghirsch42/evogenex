import pandas as pd
from scipy import stats
import numpy as np
from statsmodels.stats.multitest import multipletests


# t-test between the expression values of the treatment mice for all genes in the treatment dataset
# saves two files: 1) a file with results for all the genes in the treatment dataset and
# 2) a file with the differential expression of just the adaptive genes merged with adaptive results

ttest_results_save_file = "results/mouse_treatment/no_norm_ttest.csv"
merge_results_save_file = "results/mouse_treatment/no_norm_ttest_adaptive.csv"

print("Getting data")
expr_data = pd.read_csv("data/mouse_treatment/expr_data.csv")
response_data = pd.read_csv("data/mouse_treatment/responder_data.csv")

responders = response_data[(response_data["Response"] == "1R") | (response_data["Response"] == "2S")]
nonresponders = response_data[response_data["Response"] == "3NR"]
controls = response_data[response_data["Response"] == "4IgG"]

all_data = pd.DataFrame(columns=["regime"])
for regime in ["agg", "clade", "1_4_22", "3_10_14", "12_17_23", "5", "6_etc"]:
    evog_data = pd.read_csv("results/tpm_allrep2/orig/gene_lists/adpt_{}/all_results.csv".format(regime))
    evog_data["regime"] = regime
    all_data = all_data.append(evog_data)

print("Calculating neutral genes")
neutral_genes = all_data.groupby("gene_name").filter(lambda gene_group: not (gene_group["adaptive"] == "adaptive").any())
# print(len(neutral_genes))
b_gamma = neutral_genes[neutral_genes["brown_gamma"] < 10**10]["brown_gamma"]
# print(np.mean(np.log10(b_gamma)))
# print(np.std(np.log10(b_gamma)))
b_gamma_cutoff = np.mean(np.log10(b_gamma)) + 3*np.std(np.log10(b_gamma)) 
# print(b_gamma_cutoff)
# print(10**np.mean(np.log10(b_gamma)))
# print(10**b_gamma_cutoff)
# exit()
neutral_genes = neutral_genes[["gene_name", "brown_gamma", "qvalue"]]
neutral_genes = neutral_genes.groupby("gene_name").agg(
                                                        {
                                                            "qvalue": "min", 
                                                            "brown_gamma": "max"
                                                        }).reset_index()

neutral_genes = neutral_genes[
                    (neutral_genes["qvalue"] >= 0.5) 
                  & (np.log10(neutral_genes["brown_gamma"]) < b_gamma_cutoff)]

# print((neutral_genes))
# exit()
print("Calculating neutral expression")
neutral_expr = expr_data[expr_data["GENE"].isin(neutral_genes["gene_name"])]
# neutral_expr.to_csv("data/mouse_treatment/neutral_genes.csv", index=False)
# exit()
response_scalers = neutral_expr[responders["SYMBOL"]].mean()
nonresponse_scalers = neutral_expr[nonresponders["SYMBOL"]].mean()
control_scalers = neutral_expr[controls["SYMBOL"]].mean()
# print(response_scalers.mean())
# print(response_scalers.std())
# print(nonresponse_scalers.mean())
# print(nonresponse_scalers.std())
scalers = pd.DataFrame(response_scalers.append(nonresponse_scalers))
scalers = scalers.reset_index()
# print(scalers)
scalers.columns = ["symbol", "scale_val"]
# print(scalers["symbol"].str.split("_"))
scalers["mouse"] = scalers["symbol"].str.split("_", expand=True)[1]
# scalers.to_csv("neutral_scaler_values.csv", index=False)
# exit()

print("Running t-tests")
t_results = []

# print("expr data")
# print(expr_data)
# print("all data")
# print(all_data)
# print("evog data")
# print(evog_data)
# exit()

for gene in expr_data["GENE"].unique():
    expr = expr_data[expr_data["GENE"] == gene]
    if expr.empty:
        continue
    r_expr = expr[responders["SYMBOL"]].values[0]
    n_expr = expr[nonresponders["SYMBOL"]].values[0]
    c_expr = expr[controls["SYMBOL"]].values[0]
    # r_expr_scaled =  [r_expr[i]/response_scalers[i] for i in range(len(r_expr))]
    # n_expr_scaled =  [n_expr[i]/nonresponse_scalers[i] for i in range(len(n_expr))]
    # c_expr_scaled =  [c_expr[i]/control_scalers[i] for i in range(len(c_expr))]
    r_expr_scaled = r_expr
    n_expr_scaled = n_expr
    c_expr_scaled = c_expr
    # exit()
    ind_t = stats.ttest_ind(r_expr_scaled, n_expr_scaled)
    t_results += [[gene, ind_t[0], ind_t[1], np.mean(r_expr_scaled), np.mean(n_expr_scaled), np.mean(c_expr_scaled)]]

# print(t_results)
res_df = pd.DataFrame(t_results, columns=["gene_name", "t-statistic", "p-value", "r mean expr", "n mean expr", "c mean expr"])
na_genes = res_df[res_df["p-value"].isna()]
# print(na_genes)
res_df = res_df.dropna()
fdr_p=multipletests(pvals=res_df["p-value"], alpha=0.05, method="fdr_bh")
# print(fdr_p)
res_df["fdr p-value"] = fdr_p[1]
# print(res_df)
res_df = pd.concat([res_df, na_genes], sort=True)
res_df = res_df[["gene_name", "t-statistic", "p-value", "fdr p-value", "r mean expr", "n mean expr", "c mean expr"]]
print(res_df)
# print("Saving t-test results")
res_df.to_csv(ttest_results_save_file, index=False)
# exit()

print("Concatenating t-test data with adpativity data")
for regime in ["agg", "clade", "1_4_22", "3_10_14", "12_17_23", "5", "6_etc"]:
    evog_data = pd.read_csv("results/tpm_allrep2/orig/gene_lists/adpt_{}/all_results.csv".format(regime))
    res_df = res_df.merge(evog_data[["gene_name", 
                                        "qvalue", 
                                        "brown_gamma", 
                                        "brown_theta",
                                        "ou1_gamma", 
                                        "ou1_theta",
                                        "ou1_alpha",
                                        "ou2_gamma", 
                                        "ou2_theta",
                                        "ou2_theta_base",
                                        "ou2_alpha",
                                        "logFC"]], on="gene_name")
    res_df = res_df.rename(columns = {"qvalue": "{} q-value".format(regime),
                                    "brown_gamma": "{} brown gamma".format(regime),
                                    "brown_theta": "{} brown theta".format(regime),
                                    "ou1_gamma": "{} ou1 gamma".format(regime),
                                    "ou1_theta": "{} ou1 theta".format(regime),
                                    "ou1_alpha": "{} ou1 alpha".format(regime),
                                    "ou2_gamma": "{} ou2 gamma".format(regime),
                                    "ou2_theta": "{} ou2 theta".format(regime),
                                    "ou2_theta_base": "{} ou2 theta base".format(regime),
                                    "ou2_alpha": "{} ou2 alpha".format(regime),
                                    "logFC": "{} logFC".format(regime)})
print(res_df)
print("Saving merged results")
res_df.to_csv(merge_results_save_file, index=False)