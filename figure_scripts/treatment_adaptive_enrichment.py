import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

# Read and modify data
# df = pd.read_csv("val_results_norm.csv")
results = pd.read_csv("results/mouse_treatment/ttest_allavg_details.csv")
# df = pd.read_csv("temp2.csv")
results = results[results["fdr p-value"] < 0.05]
results["treatment de"] = ["R > NR" if g["r mean expr"] > g["n mean expr"] > 0 else "R < NR" for i, g in results.iterrows()]
# print(results)
# results = results[["gene_name", "agg q-value", "agg logFC", "clade q-value", "clade logFC",
#                    "1_4_22 q-value", "1_4_22 logFC", "3_10_14 q-value", "3_10_14 logFC","treatment_de"]]
# results = results[["gene_name", "agg q-value", "1_4_22 q-value", "3_10_14 q-value", "treatment_de"]]
# print(results)
# exit()
df = pd.DataFrame(columns = ["regime","treatment de","evogenex","count","evogenex enrichment"])
# df.columns = ["regime","treatment de","evogenex","count","evogenex enrichment"]
for regime in ["agg", "clade", "1_4_22", "3_10_14"]:
    temp = results[results["{} q-value".format(regime)].astype(float) < 0.05]
    r = 802
    nr = 414
    r_pos = len(temp[(temp["{} logFC".format(regime)] > 0) 
                   & (temp["treatment de"] == "R > NR")])
    r_neg = len(temp[(temp["{} logFC".format(regime)] < 0) 
                   & (temp["treatment de"] == "R > NR")])
    nr_pos = len(temp[(temp["{} logFC".format(regime)] > 0) 
                   & (temp["treatment de"] == "R < NR")])
    nr_neg = len(temp[(temp["{} logFC".format(regime)] < 0) 
                   & (temp["treatment de"] == "R < NR")])
    df = df.append([{"regime":regime, "treatment de": "R > NR", "evogenex":"positive", "count":r_pos, "evogenex enrichment":r_pos/r}])
    df = df.append([{"regime":regime, "treatment de": "R > NR", "evogenex":"negative", "count":r_neg, "evogenex enrichment":r_neg/r}])
    df = df.append([{"regime":regime, "treatment de": "R < NR", "evogenex":"positive", "count":nr_pos, "evogenex enrichment":nr_pos/nr}])
    df = df.append([{"regime":regime, "treatment de": "R < NR", "evogenex":"negative", "count":nr_neg, "evogenex enrichment":nr_neg/nr}])
    # df = df.append([regime, "R > NR", "negative", r_neg, r_neg/r])
    # df = df.append([regime, "R < NR", "positive", nr_pos, nr_pos/nr])
    # df = df.append([regime, "R < NR", "negative", nr_neg, nr_neg/nr])
    
    # print(temp)
    # exit()
print(df)


# temp = []
# for idx, row in results.iterrows():
#     temp.append([row["gene_name"], "agg", row["treatment de"], ])

# for regime in ["agg", "clade", "1_4_22", "3_10_14"]:
#     temp = results[(results["r mean expr"].astype(float) < results["n mean expr"].astype(float)) ]
#         # & (results["{} q-value".format(regime)] < 0.05)
#         # & (results["{} logFC".format(regime)] > 0)]
#     print(results["r mean expr"].astype(float) > results["n mean expr"].astype(float))
#     print(len(temp))

# exit()
# regime = "agg"
# df = results[results[("r mean expr").astype()]]
# df["treatment de"] = df["treatment de"].replace("positive", "R > NR")
# df["treatment de"] = df["treatment de"].replace("negative", "R < NR")
# df["evogenex"] = df["evogenex"].replace("positive", "positive adaption")
# df["evogenex"] = df["evogenex"].replace("negative", "negative adaption")
# df["regime"] = df["regime"].replace("agg", "HA-R")
# df["regime"] = df["regime"].replace("clade", "AN8")
df["regime"] = df["regime"].replace("1_4_22", "HA-S")
df["regime"] = df["regime"].replace("3_10_14", "LA-S")
# df = df[df["regime"] != "AN4"]
# print(df)
# Take only the genes that are differentially expressed and adaptive
# bp_df = df[(~df["treatment de"].isin(["none", "all"])) & (~df["evogenex"].isin(["none", "all"]))]
# print()
# exit()

# Create subplots
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, sharey=True)
bp_df = df
# Plot, separating R > N and R < N
# b1 = sns.barplot(data=bp_df[bp_df["treatment de"] == "R > NR"], x="regime", y="evogenex enrichment", hue="evogenex", ax=ax1)
# b2 = sns.barplot(data=bp_df[bp_df["treatment de"] == "R < NR"], x="regime", y="evogenex enrichment", hue="evogenex", ax=ax2)
# ax1.set_title("R > N")
# ax2.set_title("R < N")

# b1 = sns.barplot(data=bp_df[bp_df["evogenex"] == "positive"], x="treatment de", y="evogenex enrichment", hue="regime", ax=ax1)
# b1 = sns.barplot(data=bp_df[bp_df["evogenex"] == "negative"], x="treatment de", y="evogenex enrichment", hue="regime", ax=ax2)
# ax1.set_title("adaptively up-regulated")
# ax2.set_title("adaptive down-regulated")
# fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=False, sharey=True)

b1 = sns.barplot(data=bp_df[bp_df["evogenex"] == "positive"], x="regime", y="evogenex enrichment", hue="treatment de", ax=ax1)
b1 = sns.barplot(data=bp_df[bp_df["evogenex"] == "negative"], x="regime", y="evogenex enrichment", hue="treatment de", ax=ax2)
ax1.set_title("adaptively up-regulated")
ax2.set_title("adaptively down-regulated")

# b1 = sns.barplot(data=bp_df[bp_df["regime"] == "HA-R"], x="treatment de", y="evogenex enrichment", hue="evogenex", ax=ax1)
# b1 = sns.barplot(data=bp_df[bp_df["regime"] == "HA-S"], x="treatment de", y="evogenex enrichment", hue="evogenex", ax=ax2)
# b1 = sns.barplot(data=bp_df[bp_df["regime"] == "LA-S"], x="treatment de", y="evogenex enrichment", hue="evogenex", ax=ax3)
# ax1.set_title("HA-R")
# ax2.set_title("HA-S")
# ax3.set_title("LA-S")

# Create a single legend
handles = []
labels = []
axHandle, axLabel = ax1.get_legend_handles_labels()
handles.extend(axHandle)
labels.extend(axLabel)
# fig.legend(handles, labels, bbox_to_anchor=(1,.8))
# ax1.legend([],[],frameon=False)
# ax2.legend([],[],frameon=False)
# ax3.legend([],[],frameon=False)

# Create single labels
ax1.set_xlabel("")
ax1.set_ylabel("")
ax2.set_ylabel("")
# ax3.set_ylabel("")
fig.supylabel("Enrichment of adaptive genes")

# Add title
fig.suptitle("Enrichment of Adaptive Genes\nin Genes Differentially Expressed\nBetween Responders and Non-Responders")

# Save
plt.tight_layout()
# plt.savefig("figures/adaptive_enrichment.svg", format="svg")
# plt.savefig("figures/adaptive_enrichment_4.svg", format="svg")
plt.show()
