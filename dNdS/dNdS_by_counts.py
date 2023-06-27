import trisicell as tsc
import pandas as pd
import numpy as np
import csv
# from matplotlib import pyplot as plt
# import seaborn as sns

# Regime in question
regime = "all"
by_gene = False # True to give by-gene results, false to give genome-wide results

# Data paths
base_path = "results/tpm_allrep2/orig/gene_lists/adpt_{}/".format(regime)
adpt_path = base_path + "all_results.csv"
save_file = "dNdS/dnds_by_counts_results/{}.csv".format(regime)
output_arr = []

# Sublines in each regime
if regime == "har":
    sublines =  ["C18", "C15", "C11", "C16"]
# elif regime == "clade":
#     sublines = ["C13", "C18", "C15", "C11", "C16", "C8", "C20", "C7"]
elif regime == "has":
    sublines = ["C1", "C4", "C22"]
elif regime == "las":
    sublines = ["C3", "C10", "C14"]
elif regime == "all":
    sublines = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11",
                "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20",
                "C21", "C22", "C23", "C24"]
else: 
    print("Invalid regime.")
    exit()

print(regime)
print("----------")

# Load trisicell mutation data
bwes = tsc.datasets.sublines_bwes()
tsc_output = bwes.to_df(layer="trisicell_output").transpose()

# Remove clonal mutations
tsc_output = tsc_output[~(tsc_output == 1).all(axis=1)]

# Add the types of mutations
mut_kinds = tsc_output.merge(bwes.var["kind"], how="inner", left_on=tsc_output.index, right_on=bwes.var["kind"].index)
mut_kinds = mut_kinds.set_index(mut_kinds["key_0"]).drop("key_0", axis=1)

# Extract the ensemble id
mut_kinds["ensemble_id"] = mut_kinds.index.str.split(".").str[0]

# Split into synonymous and nonsynonymous
syn = mut_kinds[mut_kinds["kind"] == "synonymous SNV"]
nonsyn = mut_kinds[mut_kinds["kind"] == "nonsynonymous SNV"]

# Group by ensemble id
syn_by_gene = syn.groupby("ensemble_id").sum()
nonsyn_by_gene = nonsyn.groupby("ensemble_id").sum()
print("Total # of synonymous genes:", len(syn_by_gene))
print("Total # of nonsynonymous genes:", len(nonsyn_by_gene))
output_arr += [["Total # of synonymous genes", len(syn_by_gene)]]
output_arr += [["Total # of nonsynonymous genes", len(nonsyn_by_gene)]]
print("----------")

# Helpful variables
bkg_sublines = np.setdiff1d(syn_by_gene.columns, sublines)
n_sublines = len(sublines)
n_bkg = len(bkg_sublines)

# By gene calculations dN/dS
if by_gene:
    # Get mutation sums for the selected sublines
    sum_syn_select = syn_by_gene[sublines].sum(axis=1)
    sum_syn_select = sum_syn_select.rename("adpt").to_frame()
    sum_nonsyn_select = nonsyn_by_gene[sublines].sum(axis=1)
    sum_nonsyn_select = sum_nonsyn_select.rename("adpt").to_frame()

    # Get mutation sums for the background
    sum_syn_bkg = syn_by_gene[bkg_sublines].sum(axis=1)
    sum_syn_bkg = sum_syn_bkg.rename("bkg").to_frame()
    sum_nonsyn_bkg = nonsyn_by_gene[bkg_sublines].sum(axis=1)
    sum_nonsyn_bkg = sum_nonsyn_bkg.rename("bkg").to_frame()

    sum_syn_join = pd.merge(sum_syn_select, sum_syn_bkg, left_index=True, right_index=True)
    sum_nonsyn_join = pd.merge(sum_nonsyn_select, sum_nonsyn_bkg, left_index=True, right_index=True)

    n_sym = len(sum_syn_join)
    n_nonsym = len(sum_nonsyn_join)
    print(sum_syn_join)
    print(sum_nonsyn_join)
    exit()
    # Take only the ones that have at least one mutation on each side
    # sum_join = sum_join[(sum_join["bkg"]!=0) & (sum_join["adpt"]!=0)]
    sum_join["% adpt"] = sum_join["adpt"] / n_sublines
    sum_join["% bkg"] = sum_join["bkg"] / n_bkg
    sum_join["% bkg"] = sum_join["% bkg"].replace(0, 10**-10)
    sum_join["dN/dS"] = sum_join["% adpt"] / sum_join["% bkg"]
    results = sum_join[sum_join["dN/dS"] != 1]
    # print(len(results), "(", round(len(results)/n*100, 1), "%)")
    # print(results.index)
    # print(len(results[results["dN/dS"] > 1]), "(", round(len(results[results["dN/dS"] > 1])/n*100, 1), "%)")
    # print(len(results[results["dN/dS"] < 1]), "(", round(len(results[results["dN/dS"] < 1])/n*100, 1), "%)")
    # results.to_csv("dNdS_{}.csv".format(regime))
    # with open(save_path, "w") as f:
    #     for id in sum_join[sum_join["dN/dS"] != 1].index.to_list():
    #         f.write(id +"\n")
    print("Number with no mutated genes in selected:", len(results[results["dN/dS"] == 0]))
    print("Number dN/dS < 1:", len(results[(results["dN/dS"] < 1) & (results["dN/dS"] != 0)]))
    print("Number dN/dS > 1:", len(results[(results["dN/dS"] > 1) & (results["% bkg"] != 10**-10)]))
    print("Number with no mutated genes in background:", len(results[results["% bkg"] == 10**-10]))
    # sns.histplot(data=results, x="dN/dS")
    # plt.savefig("temp.png")s
else: # Sums over all genes and all sublines in selected or background
    # Sum selected
    sum_syn_select = syn_by_gene[syn_by_gene.columns[syn_by_gene.columns.isin(sublines)]].sum().sum()
    sum_nonsyn_select = nonsyn_by_gene[nonsyn_by_gene.columns[nonsyn_by_gene.columns.isin(sublines)]].sum().sum()
    # Sum background
    sum_syn_bkg = syn_by_gene[syn_by_gene.columns[syn_by_gene.columns.isin(bkg_sublines)]].sum().sum()
    sum_nonsyn_bkg = nonsyn_by_gene[nonsyn_by_gene.columns[nonsyn_by_gene.columns.isin(bkg_sublines)]].sum().sum()
    output_arr += [["Sum of all synonymous mutations across all genes and selected sublines", sum_syn_select]]
    output_arr += [["Sum of all synonymous mutations across all genes and background sublines", sum_syn_bkg]]
    output_arr += [["Sum of all nonsynonymous mutations across all genes and selected sublines", sum_nonsyn_select]]
    output_arr += [["Sum of all nonsynonymous mutations across all genes and background sublines", sum_nonsyn_bkg]]
    print("Sum of all synonymous mutations across all genes and selected sublines:", sum_syn_select)
    print("Sum of all synonymous mutations across all genes and background sublines:", sum_syn_bkg)
    print("Sum of all nonsynonymous mutations across all genes and selected sublines:", sum_nonsyn_select)
    print("Sum of all nonsynonymous mutations across all genes and background sublines:", sum_nonsyn_bkg)
    
    print("Nonsyn / syn for select sublines:", round(sum_nonsyn_select / sum_syn_select, 3))
    output_arr += [["Nonsyn / syn for select sublines", round(sum_nonsyn_select / sum_syn_select, 3)]]
    if n_bkg != 0:
        print("Nonsyn / syn for background sublines:", round(sum_nonsyn_bkg / sum_syn_bkg, 3))
        output_arr += [["Nonsyn / syn for background sublines", round(sum_nonsyn_bkg / sum_syn_bkg, 3)]]

    print("----------")

    # Normalize by expected value of synonymous and nonsynonymous mutations
    
    norm_sum_syn_select = sum_syn_select * 3
    norm_sum_nonsyn_select = sum_nonsyn_select * 3 / 2
    print("Normalized sum of all synonymous mutations across all genes and selected sublines:", round(norm_sum_syn_select, 3))
    print("Normalized sum of all nonsynonymous mutations across all genes and selected sublines:", round(norm_sum_nonsyn_select, 3))
    output_arr += [["Normalized sum of all synonymous mutations across all genes and selected sublines", round(norm_sum_syn_select, 3)]]
    output_arr += [["Normalized sum of all nonsynonymous mutations across all genes and selected sublines", round(norm_sum_nonsyn_select, 3)]]

    if n_bkg != 0:
        norm_sum_syn_bkg = sum_syn_bkg * 3
        norm_sum_nonsyn_bkg = sum_nonsyn_bkg * 3 / 2
        print("Normalized sum of all synonymous mutations across all genes and background sublines:", round(norm_sum_syn_bkg, 3))
        print("Normalized sum of all nonsynonymous mutations across all genes and background sublines:", round(norm_sum_nonsyn_bkg, 3))
        output_arr += [["Normalized sum of all synonymous mutations across all genes and background sublines", round(norm_sum_syn_bkg, 3)]]
        output_arr += [["Normalized sum of all nonsynonymous mutations across all genes and background sublines", round(norm_sum_nonsyn_bkg, 3)]]

    print("Normalized nonsyn / Normalized syn for select sublines:", round(norm_sum_nonsyn_select / norm_sum_syn_select, 3))
    output_arr += [["Normalized nonsyn / Normalized syn for select sublines", round(norm_sum_nonsyn_select / norm_sum_syn_select, 3)]]
    if n_bkg != 0:
        print("Normalized nonsyn / Normalized syn for background sublines:", round(norm_sum_nonsyn_bkg / norm_sum_syn_bkg, 3))    
        output_arr += [["Normalized nonsyn / Normalized syn for background sublines", round(norm_sum_nonsyn_bkg / norm_sum_syn_bkg, 3)]]

    with open(save_file, "w") as f:
        writer = csv.writer(f)
        writer.writerows(output_arr)