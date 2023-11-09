import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import math

norm_by_neut = True

df = pd.read_csv("data/mouse_treatment/Melanocytic_antigen_in_R_vs_NR.csv", header=[0,1], index_col=0)
df = df.stack().stack().reset_index()
df.columns = ["gene", "response", "mouse", "expr"]
df["response2"] = ["responder" if df.loc[i, "response"] == "R" or df.loc[i, "response"] == "SD" 
                        else "non-responder"
                        for i in range(len(df))]

gene_order = ["Ppp3ca", "Dct", "Mlana", "Pmel", "Mitf", "Rab38", "Tyr", "Tyrp1"]
df["order"] = [gene_order.index(row["gene"]) for idx, row in df.iterrows()]
df = df.sort_values(by="order")
df = df.groupby("gene", sort=False).apply(lambda x: x.sort_values(by="response2", ascending=False))
print(df)

if norm_by_neut:
    scaler_df = pd.read_csv("data/mouse_treatment/neutral_scaler_values.csv")
    scaler_df["mouse"] = scaler_df["mouse"].astype(str)
    scaled_expr = []
    for idx, row in df.iterrows():
        scale_val = scaler_df[scaler_df["mouse"] == row["mouse"]]["scale_val"].values[0]
        scaled_expr.append(row["expr"] / scale_val)
    df["expr"] = scaled_expr

print(df)
# exit()

genes = df["gene"].unique()

fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(8, 5))
for i in range(len(genes)):
    ax = axes[math.floor(i/4)][i%4]
    sns.boxplot(data=df[df["gene"] == genes[i]], x="gene", y="expr", hue="response2", ax=ax)
    ax.get_legend().remove()
    ax.set(xlabel=None)
    ax.set(ylabel=None)
    ax.set_title(genes[i])
    ax.set_xticks([])

# exit()
fig.supylabel("normalized expression")
lines_labels = axes[0][0].get_legend_handles_labels()
fig.legend(lines_labels[0], lines_labels[1],ncol=3)
plt.suptitle("Melanocytic antigen expression")
plt.tight_layout()
plt.savefig("scratch_figures/melanocytic_antigen_boxplots_code_output.svg")
plt.savefig("scratch_figures/melanocytic_antigen_boxplots_code_output.png")
plt.show()