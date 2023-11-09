import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

# Read t-test data
results = pd.read_csv("results/mouse_treatment/differential_expression/norm_ttest_adaptive.csv")
# results = pd.read_csv("results/mouse_treatment/differential_expression/no_norm_ttest_adaptive.csv")

# Take only the significantly different genes
results = results[results["fdr p-value"] < 0.05]

# Label if the responder or non-responder is larger
results["r_nr"] = ["R > NR" if g["r mean expr"] > g["n mean expr"] > 0 else "R < NR" for i, g in results.iterrows()]

# Calculate how many DEGs are higher in responder and nonresponders
r = sum(results["r_nr"] == "R > NR")
nr = sum(results["r_nr"] == "R < NR")

print(r)
print(nr)
# exit()

df = pd.DataFrame(columns = ["regime","r_nr","pos_neg_adpt","adpt_de","adpt_de/de"])

for regime in ["agg", "1_4_22", "3_10_14"]:
    # Take only adaptive results for that regime
    temp = results[results["{} q-value".format(regime)].astype(float) < 0.05]
    # Calculate the numbers
    r_pos = len(temp[(temp["{} logFC".format(regime)] > 0) 
                   & (temp["r_nr"] == "R > NR")])
    r_neg = len(temp[(temp["{} logFC".format(regime)] < 0) 
                   & (temp["r_nr"] == "R > NR")])
    nr_pos = len(temp[(temp["{} logFC".format(regime)] > 0) 
                   & (temp["r_nr"] == "R < NR")])
    nr_neg = len(temp[(temp["{} logFC".format(regime)] < 0) 
                   & (temp["r_nr"] == "R < NR")])
    df = df.append([{"regime":regime, "r_nr": "R > NR", "pos_neg_adpt":"positive", "adpt_de":r_pos, "adpt_de/de":r_pos/r}])
    df = df.append([{"regime":regime, "r_nr": "R > NR", "pos_neg_adpt":"negative", "adpt_de":r_neg, "adpt_de/de":r_neg/r}])
    df = df.append([{"regime":regime, "r_nr": "R < NR", "pos_neg_adpt":"positive", "adpt_de":nr_pos, "adpt_de/de":nr_pos/nr}])
    df = df.append([{"regime":regime, "r_nr": "R < NR", "pos_neg_adpt":"negative", "adpt_de":nr_neg, "adpt_de/de":nr_neg/nr}])
    # df = df.append([regime, "R > NR", "negative", r_neg, r_neg/r])
    # df = df.append([regime, "R < NR", "positive", nr_pos, nr_pos/nr])
    # df = df.append([regime, "R < NR", "negative", nr_neg, nr_neg/nr])
    
    # print(temp)
    # exit()
# print(df)
# exit()


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

df["regime"] = df["regime"].replace("agg", "HA-R")
df["regime"] = df["regime"].replace("1_4_22", "HA-S")
df["regime"] = df["regime"].replace("3_10_14", "LA-S")
# df = df[df["regime"] != "AN4"]
# print(df)
# Take only the genes that are differentially expressed and adaptive
# bp_df = df[(~df["treatment de"].isin(["none", "all"])) & (~df["evogenex"].isin(["none", "all"]))]
# print()
# exit()

df.to_csv("results/mouse_treatment/differential_expression/nnorm_ttest_summary.csv")
# df.to_csv("results/mouse_treatment/differential_expression/no_norm_ttest_summary.csv")
exit()

# Create subplots
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, sharey=True)
bp_df = df
bp_df["adpt_de/de"] = round(bp_df["adpt_de/de"] * 100, 1)
print(bp_df)
# Plot, separating R > N and R < N
b1 = sns.barplot(data=bp_df[bp_df["r_nr"] == "R > NR"], x="regime", y="adpt_de/de", hue="pos_neg_adpt", ax=ax1)
b2 = sns.barplot(data=bp_df[bp_df["r_nr"] == "R < NR"], x="regime", y="adpt_de/de", hue="pos_neg_adpt", ax=ax2)
ax1.set_title("R > N")
ax2.set_title("R < N")

# b1 = sns.barplot(data=bp_df[bp_df["evogenex"] == "positive"], x="treatment de", y="evogenex enrichment", hue="regime", ax=ax1)
# b1 = sns.barplot(data=bp_df[bp_df["evogenex"] == "negative"], x="treatment de", y="evogenex enrichment", hue="regime", ax=ax2)
# ax1.set_title("adaptively up-regulated")
# ax2.set_title("adaptive down-regulated")
# fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=False, sharey=True)

# b1 = sns.barplot(data=bp_df[bp_df["evogenex"] == "positive"], x="regime", y="evogenex enrichment", hue="treatment de", ax=ax1)
# b1 = sns.barplot(data=bp_df[bp_df["evogenex"] == "negative"], x="regime", y="evogenex enrichment", hue="treatment de", ax=ax2)
# ax1.set_title("adaptively up-regulated")
# ax2.set_title("adaptively down-regulated")

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

# Create single labels
ax1.set_xlabel("")
ax1.set_ylabel("")
ax2.set_ylabel("")
fig.supylabel("Enrichment of adaptive genes")

# Add title
fig.suptitle("Enrichment of Adaptive Genes\nin Genes Differentially Expressed\nBetween Responders and Non-Responders")

# Save
plt.tight_layout()
# plt.savefig("scratch_figures/adpative_in_no_norm_diff_eq_code_output.png")
# plt.savefig("scratch_figures/adpative_in_diff_eq_code_output.svg")
# plt.show()

# Save the data into a table
# bp_df.to_csv("results/mouse_treatment/adaptive_genes_in_diff_expr.csv", index=False)