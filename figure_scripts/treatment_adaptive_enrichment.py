import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

# Read and modify data
df = pd.read_csv("val_results_norm.csv")
df.columns = ["regime","treatment de","evogenex","count","evogenex enrichment"]
df["treatment de"] = df["treatment de"].replace("positive", "R > NR")
df["treatment de"] = df["treatment de"].replace("negative", "R < NR")
df["evogenex"] = df["evogenex"].replace("positive", "positive adaption")
df["evogenex"] = df["evogenex"].replace("negative", "negative adaption")
df["regime"] = df["regime"].replace("agg", "AN4")
df["regime"] = df["regime"].replace("clade", "AN8")
df["regime"] = df["regime"].replace("1_4_22", "AR")
df["regime"] = df["regime"].replace("3_10_14", "NR")

# Take only the genes that are differentially expressed and adaptive
bp_df = df[(~df["treatment de"].isin(["none", "all"])) & (~df["evogenex"].isin(["none", "all"]))]

# Create subplots
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, sharey=True)

# Plot, separating R > N and R < N
b1 = sns.barplot(data=bp_df[bp_df["treatment de"] == "R > NR"], x="regime", y="evogenex enrichment", hue="evogenex", ax=ax1)
b2 = sns.barplot(data=bp_df[bp_df["treatment de"] == "R < NR"], x="regime", y="evogenex enrichment", hue="evogenex", ax=ax2)
ax1.set_title("R > N")
ax2.set_title("R < N")

# Create a single legend
handles = []
labels = []
axHandle, axLabel = ax1.get_legend_handles_labels()
handles.extend(axHandle)
labels.extend(axLabel)
fig.legend(handles, labels, bbox_to_anchor=(1,.8))
ax1.legend([],[],frameon=False)
ax2.legend([],[],frameon=False)

# Create single labels
ax1.set_xlabel("")
ax1.set_ylabel("")
ax2.set_ylabel("")
fig.supylabel("nrichment of adaptive genes")

# Add title
fig.suptitle("Enrichment of Adaptive Genes\nin Genes Differentially Expressed\nBetween Responders and Non-Responders")

# Save
plt.tight_layout()
plt.savefig("figures/adaptive_enrichment.png")
# plt.show()
