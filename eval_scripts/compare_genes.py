import pandas as pd
import os
import itertools
import argparse

def print_counts(df):
    print("Number of genes in each list")

    gb = df.groupby(["tree", "adpt"]).count()
    print(gb["gene_name"])

def compare_trees(df):
    
    for adpt in df["adpt"].unique():
        print("Adpt: {}".format(adpt))
        tree_dfs = {}
        for tree in df["tree"].unique():
            tree_dfs[tree] = df[(df["tree"] == tree) & (df["adpt"] == adpt)]
        for t1, t2 in itertools.combinations(tree_dfs.keys(), 2):
            both = pd.merge(tree_dfs[t1], tree_dfs[t2], how="inner", on="gene_name")
            print("{}, {}: {}".format(t1, t2, len(both)))

def compare_adpt(df):

    for tree in df["tree"].unique():
        print("Tree: {}".format(tree))
        adpt_dfs = {}
        for adpt in df["adpt"].unique():
            adpt_dfs[adpt] = df[(df["tree"] == tree) & (df["adpt"] == adpt)]
        for t1, t2 in itertools.combinations(adpt_dfs.keys(), 2):
            both = pd.merge(adpt_dfs[t1], adpt_dfs[t2], how="inner", on="gene_name")
            print("{}, {}: {}".format(t1, t2, len(both)))
    
def compare_all(df):
    print("Genes in all experiments:")
    vcs = df["gene_name"].value_counts()
    print(vcs[vcs==len(df["tree"].unique())*len(df["adpt"].unique())])

def create_table(df, save_path):

    df2 = df[["ensemble_id", "gene_name"]].drop_duplicates()
    df2 = df2.reset_index(drop=True)
    df2 = df2.fillna(0)

    for tree in df["tree"].unique():
        for adpt in df["adpt"].unique():
            curr_df = df[(df["tree"] == tree) & (df["adpt"] == adpt)]
            curr_genes = curr_df["gene_name"].to_list()
            df2["{}_{}".format(tree, adpt)] = [curr_df[curr_df["gene_name"] == g]["qvalue"].values[0] if g in curr_genes else 1 for g in df2["gene_name"]]

    print(df2)

    if save_path:
        df2.to_csv(save_path, index=False)

def main(base_path):
    orig_path = base_path + "orig/gene_lists/"
    scaled_path = base_path + "lp_scaled/gene_lists/"

    df = pd.DataFrame()
    for folder in [f for f in os.listdir(orig_path) if "adpt" in f]:
        f = [f for f in os.listdir(orig_path+folder+"/") if "gene_info" in f][0]
        if "agg" in folder or "clade" in folder: sep = " "
        else: sep = ","
        temp = pd.read_csv(orig_path+folder+"/"+f, sep=sep)
        temp["tree"] = "orig"
        temp["adpt"] = folder.split("_")[-1]
        df = df.append(temp)
    for folder in [f for f in os.listdir(scaled_path) if "adpt" in f]:
        f = [f for f in os.listdir(scaled_path+folder+"/") if "gene_info" in f][0]
        if "agg" in folder or "clade" in folder: sep = " "
        else: sep = ","
        temp = pd.read_csv(scaled_path+folder+"/"+f, sep=sep)
        temp["tree"] = "scaled"
        temp["adpt"] = folder.split("_")[-1]
        df = df.append(temp)

    # print(df["qvalue"])

    print_counts(df)
    # compare_trees(df)
    # compare_adpt(df)
    # compare_all(df)
    # create_table(df, "results/tpm_allrep2/gene_comp.csv")
    # create_table(df, "temp.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_path", type=str, default="results/tpm_allrep2/")
    args = parser.parse_args()
    main(args.base_path)