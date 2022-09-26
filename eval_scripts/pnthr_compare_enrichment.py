import pandas as pd
import os
import itertools

def create_df(base_path):
    orig_path = base_path + "orig/gene_lists/"
    scaled_path = base_path + "lp_scaled/gene_lists/"

    df = pd.DataFrame()

    for folder in [f for f in os.listdir(orig_path) if "adpt" in f]:
        for t in ["bio", "cell", "mol"]:
            f = orig_path + folder + "/pnthr_" + t + ".csv"
            temp = pd.read_csv(f)
            temp["tree"] = "orig"
            temp["adpt"] = folder.split("_")[-1]
            temp["type"] = t
            df = df.append(temp)
    for folder in [f for f in os.listdir(scaled_path) if "adpt" in f]:
        for t in ["bio", "cell", "mol"]:
            f = scaled_path + folder + "/pnthr_" + t + ".csv"
            temp = pd.read_csv(f)
            temp["tree"] = "scaled"
            temp["adpt"] = folder.split("_")[-1]
            temp["type"] = t
            df = df.append(temp)
    return df

def print_counts(df):
    counts = df.groupby(["tree", "adpt", "type"]).count()["term.id"]
    print(counts)

def compare_trees(df):

    for adpt in df["adpt"].unique():
        print("Adpt: {}".format(adpt))
        tree_dfs = {}
        for tree in df["tree"].unique():
            for t in df["type"].unique():
                tree_dfs["{}_{}".format(tree,t)] = df[(df["tree"] == tree) & (df["adpt"] == adpt) & (df["type"] == t)]
        for t1, t2 in itertools.combinations(tree_dfs.keys(), 2):
            if t1.split("_")[-1] != t2.split("_")[-1]: continue # Only compare the same types
            both = pd.merge(tree_dfs[t1], tree_dfs[t2], how="inner", on="term.id")
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
            both = pd.merge(adpt_dfs[t1], adpt_dfs[t2], how="inner", on="term.id")
            print("{}, {}: {}".format(t1, t2, len(both)))
    



def compare_all(df):
    vcs = df["term.label"].value_counts()
    print(vcs[vcs==len(df["tree"].unique())*len(df["adpt"].unique())])

def create_table(df, save_file=None):
    results = pd.DataFrame(columns=["type", "GO ID", "Description", "orig_clade", "orig_agg", "scaled_clade", "scaled_agg"])
    results[["type", "GO ID", "Description"]] = df[["term.class", "term.id", "term.label"]].drop_duplicates()
    results = results.reset_index(drop=True)

    for idx, row in results.iterrows():
        oc = df[(df["term.id"] == row["GO ID"]) & (df["tree"] == "orig") & (df["adpt"] == "clade")]
        if len(oc) == 1:
            results.iloc[idx]["orig_clade"] = oc["pValue"].values[0]
        else:
            results.iloc[idx]["orig_clade"] = 1

        oa = df[(df["term.id"] == row["GO ID"]) & (df["tree"] == "orig") & (df["adpt"] == "agg")]
        if len(oa) == 1:
            results.iloc[idx]["orig_agg"] = oa["pValue"].values[0]
        else:
            results.iloc[idx]["orig_agg"] = 1

        sc = df[(df["term.id"] == row["GO ID"]) & (df["tree"] == "scaled") & (df["adpt"] == "clade")]
        if len(sc) == 1:
            results.iloc[idx]["scaled_clade"] = sc["pValue"].values[0]
        else:
            results.iloc[idx]["scaled_clade"] = 1

        sa = df[(df["term.id"] == row["GO ID"]) & (df["tree"] == "scaled") & (df["adpt"] == "agg")]
        if len(sa) == 1:
            results.iloc[idx]["scaled_agg"] = sa["pValue"].values[0]
        else:
            results.iloc[idx]["scaled_agg"] = 1
    
    print(results)

    if save_file:
        results.to_csv(save_file, index=False)



def compare_results(df):
    orig_clade = df[(df["tree"] == "orig") & (df["adpt"] == "clade")]
    orig_agg = df[(df["tree"] == "orig") & (df["adpt"] == "agg")]
    lp_clade = df[(df["tree"] == "lp") & (df["adpt"] == "clade")]
    lp_agg = df[(df["tree"] == "lp") & (df["adpt"] == "agg")]

    print("Orig clade has {} terms.".format(len(orig_clade)))
    print("Orig agg has {} terms.".format(len(orig_agg)))
    print("LP clade has {} terms.".format(len(lp_clade)))
    print("LP agg has {} terms.".format(len(lp_agg)))

    # compare_trees(orig_clade, orig_agg, lp_clade, lp_agg)
    # compare_adpt(orig_clade, orig_agg, lp_clade, lp_agg)
    compare_all(orig_clade, orig_agg, lp_clade, lp_agg)

def main():
    df = create_df("results/tpm_allrep2/")
    # print_counts(df)
    # compare_trees(df)
    # compare_adpt(df)
    compare_all(df)
    # compare_results(df)
    # df.to_csv("results/tpm_allrep2/pthr_go_term_comp.csv")
    # create_table(df, "results/tpm_allrep2/pthr_go_term_comp.csv")


if __name__ == "__main__":
    main()