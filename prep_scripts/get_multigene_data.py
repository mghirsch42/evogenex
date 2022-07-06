import trisicell as tsc
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

# TODO: Currently just taking the first 4 replicates, need to change that

wide = False
start = 0
end = 55401
inc = 1000
dtype = "tpm" # tpm or fpkm
save_path = "data/tpm_all/"

while start < end:
    curr_end = start + inc
    if curr_end > end: curr_end = end

    gene_indices = list(range(start,curr_end))
    gene_str = "{}-{}".format(start, curr_end)

    print("Running for genes " + gene_str)

    sc_data = tsc.datasets.sublines_scrnaseq()
    exp_data = sc_data["expression"]

    sc_df = pd.DataFrame()
    sc_df["subline"] = sc_data["expression"].obs["clone"]
    sc_df["subclone"] = sc_data["expression"].obs.index
    sc_df = sc_df.drop(columns=["subclone"])
    sc_df = sc_df.reset_index()


    if wide:
        exp_df = pd.DataFrame()
        for index in gene_indices:
            gene_row = {"gene": index}
            for subline in sc_df["subline"].unique():
                subclones = sc_df.loc[sc_df["subline"] == subline]["subclone"].tolist()
                subclones.sort()
                subline_row = []
                for subclone in subclones:
                    if len(subline_row) == 4:
                        break
                    rep = exp_data[exp_data.obs.index == subclone].layers[dtype]
                    if len(rep) != 1:
                        print("error", subline)
                        exit()
                    subline_row.append(rep[0][index])
                d = {"{}_R1".format(subline): subline_row[0], "{}_R2".format(subline): subline_row[1],
                    "{}_R3".format(subline): subline_row[2], "{}_R4".format(subline): subline_row[3]}
                gene_row.update(d)
            exp_df = exp_df.append(gene_row, ignore_index=True)
    else:
        exp_df = pd.DataFrame(columns=["gene", "species", "replicate", "exprval"])
        for index in gene_indices:
            for subline in sc_df["subline"].unique():
                subclones = sc_df.loc[sc_df["subline"] == subline]["subclone"].tolist()
                subclones.sort()
                rep_count = 1
                for subclone in subclones:
                    if rep_count == 4:
                        break
                    rep = exp_data[exp_data.obs.index == subclone].layers[dtype]
                    if len(rep) != 1:
                        print("error", subline)
                        exit()
                    exp_df = exp_df.append({
                        "gene": index,
                        "species": subline,
                        "replicate": rep_count,
                        "exprval": rep[0][index]
                    }, ignore_index=True)
                    rep_count+=1
    # print(exp_df)
    if gene_str == "":
        gene_str = "_".join(list(map(str, gene_indices)))
    exp_df.to_csv("{}gene_{}.csv".format(save_path, gene_str), index=False)
    # print("Saved genes " + gene_str + " to " + save_path)
    start = curr_end