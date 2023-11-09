import pandas as pd

deg = pd.read_csv("results/mouse_treatment/differential_expression/deseq2_results.csv")

deg = deg.rename(columns={"Unnamed: 0": "gene_id"})
deg[["ensemble_id", "gene_name"]] = deg["gene_id"].str.split("_", expand=True)
deg = deg.drop("gene_id", axis=1)
# print(deg)

print("Total number of genes for DE:", len(deg))
deg = deg.dropna()
print("Number of nonzero genes:", len(deg))
deg = deg[deg["padj"] < 0.05]
print("Number of DEGs:", len(deg))
print(deg)
print("Number DEGs with R > NR", len(deg[deg["log2FoldChange"] > 0]))
print("Number DEGs with R < NR", len(deg[deg["log2FoldChange"] < 0]))
print("---------")

deg_adpt = deg
for regime in ["agg", "1_4_22", "3_10_14"]:
    curr_data = pd.read_csv("results/tpm_allrep2/orig/gene_lists/adpt_{}/all_results.csv".format(regime))
    # print(regime, len(curr_data[curr_data["qvalue" < 0.05]]))

    curr_data = curr_data.drop("gene_name", axis=1)
    curr_data = curr_data.rename(columns = {"qvalue": "{} q-value".format(regime),
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
    
    deg_adpt = deg_adpt.merge(curr_data, on="ensemble_id")
deg_adpt.to_csv("results/mouse_treatment/differential_expression/deseq2_adpt_merge.csv", index=False)

# print(deg_adpt)
print("DEGs evaluated for adaptivity")
print("Total number of genes for DE that are in single cell:", len(deg_adpt))
print("Number DEGs with R > NR", len(deg_adpt[deg_adpt["log2FoldChange"] > 0]))
print("Number DEGs with R < NR", len(deg_adpt[deg_adpt["log2FoldChange"] < 0]))
print("---------")

for regime in ["agg", "1_4_22", "3_10_14"]:
    print("Regime", regime)
    curr_deg = deg_adpt[deg_adpt["{} q-value".format(regime)] < 0.05]
    print("Total DEG that are adaptive:", len(curr_deg))
    print("Number DEGs with R > NR", len(curr_deg[curr_deg["log2FoldChange"] > 0]))
    print("Number DEGs with R < NR", len(curr_deg[curr_deg["log2FoldChange"] < 0]))
    print("Number up-adaptive:", len(curr_deg[curr_deg["{} logFC".format(regime)] > 0]))
    print("Number down-adaptive:", len(curr_deg[curr_deg["{} logFC".format(regime)] < 0]))
    print("R > NR; Up-adaptive:", len(curr_deg[(curr_deg["log2FoldChange"] > 0) & (curr_deg["{} logFC".format(regime)] > 0)]))
    print("R > NR; Down-adaptive:", len(curr_deg[(curr_deg["log2FoldChange"] > 0) & (curr_deg["{} logFC".format(regime)] < 0)]))
    print("R < NR; Up-adaptive:", len(curr_deg[(curr_deg["log2FoldChange"] < 0) & (curr_deg["{} logFC".format(regime)] > 0)]))
    print("R < NR; Down-adaptive:", len(curr_deg[(curr_deg["log2FoldChange"] < 0) & (curr_deg["{} logFC".format(regime)] < 0)]))
    print("---------")
