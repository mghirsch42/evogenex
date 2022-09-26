import pandas as pd

def create_df(base_path):
    orig_path = base_path + "orig/gene_lists/kegg/"
    scaled_path = base_path + "lp_scaled/gene_lists/kegg/"

    df = pd.DataFrame()

    orig_clade = pd.read_csv(orig_path+"adpt_clade.csv")
    orig_clade[["tree", "adpt"]] = ["orig", "clade"]

    orig_agg = pd.read_csv(orig_path+"adpt_agg.csv")
    orig_agg[["tree", "adpt"]] = ["orig", "agg"]

    scaled_clade = pd.read_csv(scaled_path+"adpt_clade.csv")
    scaled_clade[["tree", "adpt"]] = ["scaled", "clade"]

    scaled_agg = pd.read_csv(scaled_path+"adpt_agg.csv")
    scaled_agg[["tree", "adpt"]] = ["scaled", "agg"]

    df = df.append([orig_clade, orig_agg, scaled_clade, scaled_agg])

    # print(df)
    return df

def print_counts(df):
    counts = df.groupby(["tree", "adpt"]).count()["ID"]
    print(counts)

def compare_trees(orig_clade, orig_agg, lp_clade, lp_agg):

    # Orig clade vs LP clade
    print("Orig clade vs LP clade")
    both = pd.merge(orig_clade, lp_clade, how="inner", on="ID")
    # print(both["term.label"])
    print("Shared:", len(both))

    orig_only = pd.merge(orig_clade, both, how="left", on="ID", indicator=True)
    orig_only = orig_only[orig_only["_merge"] == "left_only"]
    # print(orig_only["term.label"])
    print("Orig only:", len(orig_only))

    # Orig agg vs LP agg
    print("Orig agg vs LP agg")
    both = pd.merge(orig_agg, lp_agg, how="inner", on="ID")
    print("Shared:", len(both))
    # print(both["term.label"])

    orig_only = pd.merge(orig_agg, both, how="left", on="ID", indicator=True)
    orig_only = orig_only[orig_only["_merge"] == "left_only"]
    # print(orig_only["term.label"])
    print("Orig only:", len(orig_only))

def compare_adpt(orig_clade, orig_agg, lp_clade, lp_agg):

    # Orig clade vs LP clade
    print("Orig clade vs Orig agg")
    both = pd.merge(orig_clade, orig_agg, how="inner", on="ID")
    # print(both["term.label"])
    print("Shared:", len(both))

    clade_only = pd.merge(orig_clade, both, how="left", on="ID", indicator=True)
    clade_only = clade_only[clade_only["_merge"] == "left_only"]
    # print(orig_only["term.label"])
    print("Clade only:", len(clade_only))

    # Orig agg vs LP agg
    print("LP clade vs LP agg")
    both = pd.merge(lp_clade, lp_agg, how="inner", on="ID")
    print("Shared:", len(both))
    # print(both["term.label"])

    clade_only = pd.merge(orig_agg, both, how="left", on="ID", indicator=True)
    clade_only = clade_only[clade_only["_merge"] == "left_only"]
    # print(orig_only["term.label"])
    print("Clade only:", len(clade_only))

def compare_all(orig_clade, orig_agg, lp_clade, lp_agg):
    # both = pd.merge(orig_clade, orig_agg, how="inner", on="term.label", suffixes=["", "_y"])
    # print(len(both))
    # both = pd.merge(both, lp_clade, how="inner", on="term.label", suffixes=["", "_y"])
    # print(len(both))
    # both = pd.merge(both, lp_agg, how="inner", on="term.label", suffixes=["", "_y"])
    # print(len(both))
    # print(both["term.id"])
    # print(both["term.label"])
    # print(both["term.class"])

    all = pd.concat([orig_clade, orig_agg, lp_clade, lp_agg])
    vcs = all["Description"].value_counts()
    print(vcs[vcs==4])

def create_table(df, save_file=None):
    results = pd.DataFrame(columns=["GO ID", "Description", "orig_clade", "orig_agg", "scaled_clade", "scaled_agg", "genes"])
    results[["GO ID", "Description"]] = df[["ID", "Description"]].drop_duplicates()
    results = results.reset_index(drop=True)

    for idx, row in results.iterrows():
        oc = df[(df["ID"] == row["GO ID"]) & (df["tree"] == "orig") & (df["adpt"] == "clade")]
        if len(oc) == 1:
            results.iloc[idx]["orig_clade"] = oc["qvalue"].values[0]
        else:
            results.iloc[idx]["orig_clade"] = 1

        oa = df[(df["ID"] == row["GO ID"]) & (df["tree"] == "orig") & (df["adpt"] == "agg")]
        if len(oa) == 1:
            results.iloc[idx]["orig_agg"] = oa["qvalue"].values[0]
        else:
            results.iloc[idx]["orig_agg"] = 1

        sc = df[(df["ID"] == row["GO ID"]) & (df["tree"] == "scaled") & (df["adpt"] == "clade")]
        if len(sc) == 1:
            results.iloc[idx]["scaled_clade"] = sc["qvalue"].values[0]
        else:
            results.iloc[idx]["scaled_clade"] = 1

        sa = df[(df["ID"] == row["GO ID"]) & (df["tree"] == "scaled") & (df["adpt"] == "agg")]
        if len(sa) == 1:
            results.iloc[idx]["scaled_agg"] = sa["qvalue"].values[0]
        else:
            results.iloc[idx]["scaled_agg"] = 1
        results.iloc[idx]["genes"] = df[df["ID"] == row["GO ID"]]["geneID"].values[0]
    
    print(results)

    if save_file:
        results.to_csv(save_file, index=False)

def compare_results(df):
    orig_clade = df[(df["tree"] == "orig") & (df["adpt"] == "clade")]
    orig_agg = df[(df["tree"] == "orig") & (df["adpt"] == "agg")]
    scaled_clade = df[(df["tree"] == "scaled") & (df["adpt"] == "clade")]
    scaled_agg = df[(df["tree"] == "scaled") & (df["adpt"] == "agg")]

    print("Orig clade has {} terms.".format(len(orig_clade)))
    print("Orig agg has {} terms.".format(len(orig_agg)))
    print("Scaled clade has {} terms.".format(len(scaled_clade)))
    print("Scaled agg has {} terms.".format(len(scaled_agg)))
    print("==========")
    compare_trees(orig_clade, orig_agg, scaled_clade, scaled_agg)
    print("==========")
    compare_adpt(orig_clade, orig_agg, scaled_clade, scaled_agg)
    compare_all(orig_clade, orig_agg, scaled_clade, scaled_agg)

def main():
    df = create_df("results/tpm_allrep2/")
    # print(df.columns)
    # print_counts(df)
    # compare_results(df)
    # create_table(df, None)
    create_table(df, "results/tpm_allrep2/cp_kegg_term_comp.csv")


if __name__ == "__main__":
    main()