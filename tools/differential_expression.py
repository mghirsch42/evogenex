import pandas as pd
import os
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import numpy as np

base_path = "data/tpm_allrep2_raw/"
regime ="c_3_10_14"
regime_file = "regime_files/resolved/{}.csv".format(regime)
save_path = "tables/tpm_allrep2/dge2/{}.csv".format(regime)

regime_df = pd.read_csv(regime_file)
# print(regime_df)

r1_taxa = regime_df[(regime_df["regime"] == "nonaggressive") &  (regime_df["node2"].isna())]["node"].to_list()
r2_taxa = regime_df[(regime_df["regime"] == "aggressive") &  (regime_df["node2"].isna())]["node"].to_list()

# print(r1_taxa)
# print(r2_taxa)

df = pd.DataFrame()

for f in os.listdir(base_path):
    temp = pd.read_csv(base_path + f)
    df = df.append(temp, ignore_index = True)

df["aggressive"] = df["species"].isin(r2_taxa)

# print(df)
r1_data = df[df["species"].isin(r1_taxa)]
# print(r1_data)
r2_data = df[df["species"].isin(r2_taxa)]
# print(r2_data)

# ind_t = stats.ttest_ind(r1_data["exprval"].to_list(), r2_data["exprval"].to_list())
# print(ind_t)
# # sns.boxplot(data = df, x="aggressive", y="exprval")
# # plt.show()

print(round(r1_data["exprval"].mean(), 2))
print(round(r2_data["exprval"].mean(), 2))
# print(round(r1_data["exprval"].std(), 4))
# print(round(r2_data["exprval"].std(), 4))
exit()
results = []

for gene in df["gene"].unique():
    ensemble_id, gene_name = gene.split("_")
    curr_data = df[df["gene"] == gene]
    r1_data = curr_data[curr_data["species"].isin(r1_taxa)]
    r2_data = curr_data[curr_data["species"].isin(r2_taxa)]
    
    r1_mean = r1_data["exprval"].mean()
    r2_mean = r2_data["exprval"].mean()
    r1_std = r1_data["exprval"].std()
    r2_std = r2_data["exprval"].std()

    ind_t = stats.ttest_ind(r1_data["exprval"].to_list(), r2_data["exprval"].to_list())
    logFC = np.log2(r2_mean/r1_mean)
    results.append([ensemble_id, gene_name, r1_mean, r2_mean, r2_mean-r1_mean, r1_std, r2_std, logFC, ind_t[1]])


results_df = pd.DataFrame(data=results, columns=["ensemble_id", "gene_name", "base_mean", "adpt_mean", "mean_diff(adpt-base)", "base_std", "adpt_std", "logFC", "tt_pvalue"])
fdr_p=multipletests(pvals=results_df["tt_pvalue"], alpha=0.05, method="fdr_bh")
results_df["fdr_p"] = fdr_p[1]
print(results_df)
results_df.to_csv(save_path, index=False)
