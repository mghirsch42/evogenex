import trisicell as tsc
import pandas as pd
import numpy as np

# Regime in question
regime = "1_4_22"
by_gene = False # True to give by-gene results, false to give summary results

# Data paths
base_path = "results/tpm_allrep2/orig/gene_lists/adpt_{}/".format(regime)
adpt_path = base_path + "all_results.csv"
save_path = "{}_{}.csv".format("mutation_info", regime)

# Sublines in each regime
if regime == "agg":
    sublines =  ["C18", "C15", "C11", "C16"]
elif regime == "clade":
    sublines = ["C13", "C18", "C15", "C11", "C16", "C8", "C20", "C7"]
elif regime == "1_4_22":
    sublines = ["C1", "C4", "C22"]
elif regime == "3_10_14":
    sublines = ["C3", "C10", "C14"]
else: 
    print("Invalid regime.")
    exit()

# Load trisicell mutation data
bwes = tsc.datasets.sublines_bwes()
tsc_output = bwes.to_df(layer="trisicell_output").transpose()
tsc_output["ensemble_id"] = tsc_output.index.str.split(".").str[0]
tsc_by_gene = tsc_output.groupby("ensemble_id").sum()
print("# of mutated genes", len(tsc_by_gene))

bkg_sublines = np.setdiff1d(tsc_by_gene.columns, sublines)
n_sublines = len(sublines)
n_bkg = len(bkg_sublines)

if by_gene:
    # By gene calculations dN/dS
    sum_sublines = tsc_by_gene[sublines].sum(axis=1)
    sum_sublines = sum_sublines.rename("adpt").to_frame()
    sum_bkg = tsc_by_gene[bkg_sublines].sum(axis=1)
    sum_bkg = sum_bkg.rename("bkg").to_frame()
    sum_join = pd.merge(sum_sublines, sum_bkg, left_index=True, right_index=True)
    sum_join = sum_join[(sum_join["bkg"]!=0) & (sum_join["adpt"]!=0)]
    sum_join["% adpt"] = sum_join["adpt"] / n_sublines
    sum_join["% bkg"] = sum_join["bkg"] / n_bkg
    sum_join["dN/dS"] = sum_join["% adpt"] / sum_join["% bkg"]
    print(len(sum_join[sum_join["dN/dS"] != 1]))
    print(sum_join[sum_join["dN/dS"] != 1].index)
    with open(save_path, "w") as f:
        for id in sum_join[sum_join["dN/dS"] != 1].index.to_list():
            f.write(id +"\n")
else:
    # Sums over all genes
    sum_sublines = tsc_by_gene[tsc_by_gene.columns[tsc_by_gene.columns.isin(sublines)]].sum().sum()
    sum_bkg = tsc_by_gene[tsc_by_gene.columns[tsc_by_gene.columns.isin(bkg_sublines)]].sum().sum()
    norm_sublines = sum_sublines / n_sublines
    norm_bkg = sum_bkg / n_bkg

    print("Sublines ({}) sums across all genes:".format(n_sublines))
    print("{}, {}".format(sum_sublines, round(norm_sublines,2)))
    print("Background ({}) sums across all genes:".format(n_bkg))
    print("{}, {}".format(sum_bkg, round(norm_bkg, 2)))
    print("norm sublines / norm background:")
    print(round(norm_sublines / norm_bkg, 3))
