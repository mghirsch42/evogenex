import trisicell as tsc
import pandas as pd
import numpy as np

# Just count # of adaptive genes with mutations

# regime = "agg"
# regime = "1_4_22"
regime = "3_10_14"

# How to count mutations in each group
# any - count genes where any sublines in that group are mutated
# all - count genes where all sublines in that group are mutated
# none - count genes where no sublines in that group are mutated
adpt_criteria = "none"
bkg_criteria = "any"

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

# Load mutation data
bwes = tsc.datasets.sublines_bwes()
tsc_output = bwes.to_df(layer="trisicell_output").transpose()
# Remove clonal mutations
tsc_output = tsc_output[~(tsc_output == 1).all(axis=1)]
# Label by gene
tsc_output["ensemble_id"] = tsc_output.index.str.split(".").str[0]
# Group by gene
tsc_by_gene = tsc_output.groupby("ensemble_id").sum()
# Make mutations binary per gene rather than sum
tsc_by_gene[tsc_by_gene > 0] = 1
print("# of mutated genes:", len(tsc_by_gene))

# Load adaptive genes for selected regime
adpt_res = pd.read_csv(adpt_path)
adpt_res = adpt_res[adpt_res["adaptive"] == "adaptive"]
adpt_res["ensemble_id"] = adpt_res["ensemble_id"].str.split(".").str[0]
adpt_up = adpt_res[adpt_res["logFC"] > 0]["ensemble_id"]
adpt_down = adpt_res[adpt_res["logFC"] < 0]["ensemble_id"]
print("# of adaptive genes:", len(adpt_res))
print("# of upwardly adaptive genes:", len(adpt_up))
print("# of downwardly adaptive genes:", len(adpt_down))
print("----------")

# Get background sublines and lengths
bkg_sublines = np.setdiff1d(tsc_by_gene.columns, sublines)
n_sublines = len(sublines)
n_bkg = len(bkg_sublines)

# Get genes that are adative and mutated
tsc_adpt = tsc_by_gene[tsc_by_gene.index.isin(adpt_res["ensemble_id"])]
print("# genes both adaptive and mutated:", len(tsc_adpt))
print("# genes both upwardly adaptive and mutated:", len(tsc_adpt[tsc_adpt.index.isin(adpt_up)]))
print("# genes both downwardly adaptive and mutated:", len(tsc_adpt[tsc_adpt.index.isin(adpt_down)]))
print("----------")

# Separate out the sublines
tsc_sublines = tsc_adpt[sublines]
tsc_bkg = tsc_adpt[bkg_sublines]
if adpt_criteria == "any":
    select_sublines = tsc_sublines[(tsc_sublines!= 0).any(axis=1)]
if adpt_criteria == "all":
    select_sublines = tsc_sublines[(tsc_sublines!= 0).all(axis=1)]
if adpt_criteria == "none":
    select_sublines = tsc_sublines[(tsc_sublines == 0).all(axis=1)]
if bkg_criteria == "any":
    select_bkg = tsc_bkg[(tsc_bkg != 0).any(axis=1)]
if bkg_criteria == "all":
    select_bkg = tsc_bkg[(tsc_bkg != 0).all(axis=1)]
if bkg_criteria == "none":
    select_bkg = tsc_bkg[(tsc_bkg == 0).all(axis=1)]

print("# selected genes in sublines:", len(select_sublines))
print("# selected upwardly adaptive genes in sublines:", len(select_sublines[select_sublines.index.isin(adpt_up)]))
print("# selected downwardly adaptive genes in sublines:", len(select_sublines[select_sublines.index.isin(adpt_down)]))
print("----------")

print("# selected genes in background:", len(select_bkg))
print("# selected upwardly adaptive genes in background:", len(select_bkg[select_bkg.index.isin(adpt_up)]))
print("# selected downwardly adaptive genes in background:", len(select_bkg[select_bkg.index.isin(adpt_down)]))
print("----------")

merged = pd.merge(select_sublines, select_bkg, left_on=select_sublines.index, right_on=select_bkg.index)
print("# genes in intersection of selections:", len(merged))
print("# upwardly adaptive genes in intersection of selections:", len(merged[merged["key_0"].isin(adpt_up)]))
print("# upwardly adaptive genes in intersection of selections:", len(merged[merged["key_0"].isin(adpt_down)]))
