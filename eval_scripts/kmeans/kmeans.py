import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
import seaborn as sns
from sklearn import cluster
from scipy.cluster.hierarchy import linkage, fclusterdata
from sklearn.metrics import silhouette_score

# Only adaptive genes in these files
results_dir = "results/tpm_allrep2/orig/gene_lists/adpt_{}/gene_info.csv"
regimes = ["agg", "1_4_22", "3_10_14"]

df = pd.DataFrame(columns=["gene_name", "regime", "v"])

for regime in regimes:
    temp = pd.read_csv(results_dir.format(regime))
    print(regime, len(temp))
    temp["regime"] = regime
    temp["v"] = np.where(temp["logFC"] < 0, np.log2(temp["qvalue"]), -1 * np.log2(temp["qvalue"])) # all q values will be < 1, so log will be negative, we want negative when log fold is negative
    df = df.append(temp[["gene_name", "regime", "v"]])
print(len(df))
# exit()

df2 = df.pivot(index = "gene_name", columns="regime", values="v")
df2 = df2.fillna(0) # genes that aren't adaptive in that regime are assigned 0
print(len(df2))
exit()
# sil = []
# kmin = 2
# kmax = 20

# for k in range(kmin, kmax):
#     preds = cluster.KMeans(n_clusters=k, random_state=0).fit(df2[["agg", "1_4_22", "3_10_14"]])
#     sil.append(silhouette_score(df2[["agg", "1_4_22", "3_10_14"]], preds.labels_, metric = "euclidean"))

# print(sil)
# print(np.argmax(sil))
# plt.plot(list(range(kmin, kmax)), sil)
# plt.show()

# best silouette was 12
preds = cluster.KMeans(n_clusters=12, random_state=0).fit_predict(df2[["agg", "1_4_22", "3_10_14"]])
df2["clusters"] = preds
df2 = df2.sort_values("clusters")
df2.to_csv("kmeans12_qval_clusts.csv")
exit()
clusts = pd.DataFrame(df2["clusters"])

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
sns.heatmap(data = df2, ax=ax1, cmap="seismic",  cbar_kws = {"location":"left"})
sns.heatmap(data = clusts, ax=ax2, cmap="tab20", cbar_kws = {"ticks": np.linspace(0, 11, 12)})
plt.yticks([])

# plt.savefig("kmeans12_qval.svg")
plt.show()