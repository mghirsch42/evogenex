import trisicell as tsc
import pandas as pd
import numpy as np

# TODO: Currently just taking the first 4 replicates, need to change that

gene_col = 6812
gene_name = ""
dtype = "tpm" # tpm or fpkm
save_path = "data/tpm_all/"


sc_data = tsc.datasets.sublines_scrnaseq()
exp_data = sc_data["expression"]
mut_data = sc_data["mutation"]

if gene_name == "" and gene_col != -1:
    gene_name = exp_data.var.index[gene_col]
if gene_name != "" and gene_col == -1:
    gene_col = np.where(exp_data.var.index == gene_name)
else:
    print("Give either gene name or gene column, not both")


sc_df = pd.DataFrame()
sc_df["subline"] = sc_data["expression"].obs["clone"]
sc_df["subclone"] = sc_data["expression"].obs.index
sc_df = sc_df.drop(columns=["subclone"])
sc_df = sc_df.reset_index()

exp_df = pd.DataFrame(columns=["species", "R1", "R2", "R3", "R4"])

for subline in sc_df["subline"].unique():
    subclones = sc_df.loc[sc_df["subline"] == subline]["subclone"].tolist()
    subclones.sort()
    row = []
    for subclone in subclones:
        if len(row) == 4:
            break
        rep = exp_data[exp_data.obs.index == subclone].layers[dtype]
        if len(rep) != 1:
            print("error", subline)
            exit()
        row.append(rep[0][gene_col])
    d = {"species": subline, "R1": row[0], "R2": row[1], "R3": row[2], "R4": row[3]}
    exp_df = exp_df.append(d, ignore_index=True)

exp_df.to_csv("{}gene_{}.csv".format(save_path, gene_col), index=False)
