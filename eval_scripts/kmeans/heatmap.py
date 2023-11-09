import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import cm
import seaborn as sns
from sklearn import cluster
from scipy.cluster.hierarchy import linkage, fclusterdata

results_dir = "results/tpm_allrep2/orig/gene_lists/adpt_{}/gene_info.csv"
regimes = ["agg", "1_4_22", "3_10_14"]

# df = pd.DataFrame(columns=["gene_name", "regime", "v"])

# for regime in regimes:
#     temp = pd.read_csv(results_dir.format(regime))
#     temp["regime"] = regime
#     temp["v"] = np.where(temp["logFC"] < 0, np.log2(temp["qvalue"]), -1 * np.log2(temp["qvalue"])) # all q values will be < 1, so log will be negative, we want negative when log fold is negative
#     df = df.append(temp[["gene_name", "regime", "v"]])

# df2 = df.pivot(index = "gene_name", columns="regime", values="v")
# df2 = df2.fillna(0)
# preds = cluster.KMeans(n_clusters=10).fit_predict(df2[["agg", "1_4_22", "3_10_14"]])
# df2["clusters"] = preds
df2 = pd.read_csv("kmeans12_qval_clusts_renamed.csv")
df2 = df2.set_index("gene_name", drop=True)
df2 = df2[regimes + ["clusters"]]  # this will order
df2 = df2.sort_values("clusters")

# order clusters as we want
# df2["order"] = ""
# order = [2, 4, 5, 7, 0, 11, 1, 6, 8, 10, 9, 3]
# x = np.linspace(0, 11, 12)
# order_pairs = list(zip(order, x))
# for i in range(len(order_pairs)):
#     df2.loc[df2["clusters"] == order_pairs[i][0],"order"] = order_pairs[i][1]
# df2 = df2.sort_values("order", axis=0)

print(df2)
# exit()

clusts = pd.DataFrame(df2["clusters"])
# df2.to_csv("kmeans10_qval_clusts.csv")
# exit()
# df2 = df2.drop(["clusters", "order"], axis=1)
df2 = df2.drop(["clusters"], axis=1)
# print((clusts))
# exit()
# print(df2)
cmap = cm.get_cmap("tab20", 12)

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
sns.heatmap(data = df2, ax=ax1, cmap="seismic",  cbar_kws = {"location":"left"})
sns.heatmap(data = clusts, ax=ax2, cmap=cmap)
# print(ax2.gci().colorbar)
# exit()
plt.yticks([])
plt.savefig("kmeans12_qval_renamed2.svg")
plt.show()

exit()

clusts = fclusterdata(df2, t=1.1547)
print(clusts)
print(len(np.unique(clusts)))

ns = np.linspace(0,len(mcolors.CSS4_COLORS)-1,len(np.unique(clusts))+1,dtype=int)
print(ns)
color_list = np.array(list(mcolors.CSS4_COLORS.keys()))[ns]
# print(color_list)

# print(color_list)
# lut = dict(zip(np.unique(clusts), "rbg"))
# print(lut)
# row_colors = clusts.map(lut)
# print(row_colors)
# sns.clustermap(iris, row_colors=row_colors)

# exit()
clust_df = pd.DataFrame(index=df2.index, data=clusts, columns=["cluster"])
clust_colors = color_list[clust_df["cluster"]]
print(clust_df)
print(clust_colors)

g = sns.clustermap(df2, row_colors=clust_colors)
print(g.dendrogram_col.linkage)
plt.show()