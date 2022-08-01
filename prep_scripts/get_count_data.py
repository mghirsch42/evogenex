import trisicell as tsc
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

# TODO: Currently just taking the first 4 replicates, need to change that

dtype = "tpm" # tpm or fpkm
save_file = "data/{}_count/{}_count.csv".format(dtype, dtype)


sc_data = tsc.datasets.sublines_scrnaseq()
exp_data = sc_data["expression"]

sc_df = pd.DataFrame()
sc_df["subline"] = sc_data["expression"].obs["clone"]
sc_df["subclone"] = sc_data["expression"].obs.index
sc_df = sc_df.drop(columns=["subclone"])
sc_df = sc_df.reset_index()

exp_df = pd.DataFrame(columns=["gene", "species", "replicate", "exprval"])
for subline in sc_df["subline"].unique():
    subclones = sc_df.loc[sc_df["subline"] == subline]["subclone"].tolist()
    subclones.sort()
    rep_count = 1
    for subclone in subclones:
        if rep_count > 4:
            break
        rep = exp_data[exp_data.obs.index == subclone].layers[dtype]
        if len(rep) != 1:
            print("error", subline)
            exit()
        exp_df = exp_df.append({
            "gene": 0,
            "species": subline,
            "replicate": rep_count,
            "exprval": np.count_nonzero(rep[0])
        }, ignore_index=True)
        rep_count+=1

exp_df.to_csv(save_file, index=False)
