import pandas as pd
import os
import itertools

def create_df(base_path):
    orig_path = base_path + "orig/gene_lists/"
    scaled_path = base_path + "lp_scaled/gene_lists/"

    df = pd.DataFrame()

    for folder in [f for f in os.listdir(orig_path) if "adpt" in f]:
        f = orig_path + folder + "/" + "cp_go.csv"
        temp = pd.read_csv(f)
        temp["tree"] = "orig"
        temp["adpt"] = folder.split("_")[-1]
        df = df.append(temp)
    for folder in [f for f in os.listdir(scaled_path) if "adpt" in f]:
        f = scaled_path + folder + "/" + "cp_go.csv"
        temp = pd.read_csv(f)
        temp["tree"] = "scaled"
        temp["adpt"] = folder.split("_")[-1]
        df = df.append(temp)
    # print(df)
    return df

def print_counts(df):
    counts = df.groupby(["tree", "adpt", "type"]).count()["ID"]
    print(counts)

def compare_trees(df):

    for adpt in df["adpt"].unique():
        print("Adpt: {}".format(adpt))
        tree_dfs = {}
        for tree in df["tree"].unique():
            for t in df["type"].unique():
                tree_dfs["{}_{}".format(tree,t)] = df[(df["tree"] == tree) & (df["adpt"] == adpt) & (df["type"] == t)]
        print(tree_dfs)
        for t1, t2 in itertools.combinations(tree_dfs.keys(), 2):
            if t1.split("_")[-1] != t2.split("_")[-1]: continue # Only compare the same types
            both = pd.merge(tree_dfs[t1], tree_dfs[t2], how="inner", on="ID")
            print("{}, {}: {}".format(t1, t2, len(both)))


def compare_adpt(df):
    for tree in df["tree"].unique():
        print("Tree: {}".format(tree))
        adpt_dfs = {}
        for adpt in df["adpt"].unique():
            for t in df["type"].unique():
                adpt_dfs["{}_{}".format(adpt, t)] = df[(df["tree"] == tree) & (df["adpt"] == adpt) & (df["type"] == t)]
        for t1, t2 in itertools.combinations(adpt_dfs.keys(), 2):
            if t1.split("_")[-1] != t2.split("_")[-1]: continue # Only compare the same types
            both = pd.merge(adpt_dfs[t1], adpt_dfs[t2], how="inner", on="ID")
            print("{}, {}: {}".format(t1, t2, len(both)))
    

def compare_all(df):
    vcs = df["Description"].value_counts()
    print(vcs[vcs==len(df["tree"].unique())*len(df["adpt"].unique())])

def create_table(df, save_file=None):
    exps = []
    for tree in df["tree"].unique():
        for adpt in df["adpt"].unique():
            exps.append("{}_{}".format(tree, adpt))
    results = pd.DataFrame(columns=["type", "GO ID", "Description"] + exps + ["genes"])
    results[["type", "GO ID", "Description", "genes"]] = df[["type", "ID", "Description", "geneID"]].drop_duplicates()
    results = results.reset_index(drop=True)

    for idx, row in results.iterrows():
        for tree in df["tree"].unique():
            for adpt in df["adpt"].unique():
                terms = df[(df["ID"] == row["GO ID"]) & (df["tree"] == tree) & (df["adpt"] == adpt)]
                if len(terms) == 1:
                    results.iloc[idx]["{}_{}".format(tree, adpt)] = terms["qvalue"].values[0]
                else:
                    results.iloc[idx]["{}_{}".format(tree, adpt)] = 1

    print(results)

    if save_file:
        results.to_csv(save_file, index=False)


def main():
    df = create_df("results/tpm_allrep2/")
    # print_counts(df)
    # compare_trees(df)
    # compare_adpt(df)
    # compare_all(df)
    # create_table(df, None)
    create_table(df, "results/tpm_allrep2/cp_go_term_comp.csv")


if __name__ == "__main__":
    main()