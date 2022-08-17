import pandas as pd

def create_df(base_path):
    orig_path = base_path + "orig/gene_lists/"
    lp_path = base_path + "lp_scaled/gene_lists/"

    df = pd.DataFrame()

    orig_clade_bio = pd.read_csv(orig_path+"adpt_clade_bio.csv")
    orig_clade_bio[["tree", "adpt", "term.class"]] = ["orig", "clade", "bio"]
    orig_clade_cell = pd.read_csv(orig_path+"adpt_clade_cell.csv")
    orig_clade_cell[["tree", "adpt", "term.class"]] = ["orig", "clade", "mol"]
    orig_clade_mol = pd.read_csv(orig_path+"adpt_clade_mol.csv")
    orig_clade_mol[["tree", "adpt", "term.class"]] = ["orig", "clade", "mol"]

    orig_agg_bio = pd.read_csv(orig_path+"adpt_agg_bio.csv")
    orig_agg_bio[["tree", "adpt", "term.class"]] = ["orig", "agg", "bio"]
    orig_agg_cell = pd.read_csv(orig_path+"adpt_agg_cell.csv")
    orig_agg_cell[["tree", "adpt", "term.class"]] = ["orig", "agg", "cell"]
    orig_agg_mol = pd.read_csv(orig_path+"adpt_agg_mol.csv")
    orig_agg_mol[["tree", "adpt", "term.class"]] = ["orig", "agg", "mol"]

    lp_clade_bio = pd.read_csv(lp_path+"adpt_clade_bio.csv")
    lp_clade_bio[["tree", "adpt", "term.class"]] = ["lp", "clade", "bio"]
    lp_clade_cell = pd.read_csv(lp_path+"adpt_clade_cell.csv")
    lp_clade_cell[["tree", "adpt", "term.class"]] = ["lp", "clade", "cell"]
    lp_clade_mol = pd.read_csv(lp_path+"adpt_clade_mol.csv")
    lp_clade_mol[["tree", "adpt", "term.class"]] = ["lp", "clade", "mol"]
    
    lp_agg_bio = pd.read_csv(lp_path+"adpt_agg_bio.csv")
    lp_agg_bio[["tree", "adpt", "term.class"]] = ["lp", "agg", "bio"]
    lp_agg_cell = pd.read_csv(lp_path+"adpt_agg_cell.csv")
    lp_agg_cell[["tree", "adpt", "term.class"]] = ["lp", "agg", "cell"]
    lp_agg_mol = pd.read_csv(lp_path+"adpt_agg_mol.csv")
    lp_agg_mol[["tree", "adpt", "term.class"]] = ["lp", "agg", "mol"]

    df = df.append([orig_clade_bio, orig_clade_cell, orig_clade_mol, 
                    orig_agg_bio, orig_agg_cell, orig_agg_mol,
                    lp_clade_bio, lp_clade_cell, lp_clade_mol,
                    lp_agg_bio, lp_agg_cell, lp_agg_mol])

    # print(df)
    return df

def print_counts(df):
    counts = df.groupby(["tree", "adpt", "term.class"]).count()["term.label"]
    print(counts)

def compare_trees(orig_clade, orig_agg, lp_clade, lp_agg):

    # Orig clade vs LP clade
    print("Orig clade vs LP clade")
    both = pd.merge(orig_clade, lp_clade, how="inner", on="term.label")
    # print(both["term.label"])
    print(len(both))

    orig_only = pd.merge(orig_clade, both, how="left", on="term.label", indicator=True)
    orig_only = orig_only[orig_only["_merge"] == "left_only"]
    # print(orig_only["term.label"])

    # Orig agg vs LP agg
    print("Orig agg vs LP agg")
    both = pd.merge(orig_agg, lp_agg, how="inner", on="term.label")
    print(len(both))
    # print(both["term.label"])

    orig_only = pd.merge(orig_agg, both, how="left", on="term.label", indicator=True)
    orig_only = orig_only[orig_only["_merge"] == "left_only"]
    # print(orig_only["term.label"])
    # print(len(orig_only))

def compare_adpt(orig_clade, orig_agg, lp_clade, lp_agg):

    # Orig clade vs LP clade
    print("Orig clade vs Orig agg")
    both = pd.merge(orig_clade, orig_agg, how="inner", on="term.label")
    # print(both["term.label"])
    print(len(both))

    clade_only = pd.merge(orig_clade, both, how="left", on="term.label", indicator=True)
    clade_only = clade_only[clade_only["_merge"] == "left_only"]
    # print(orig_only["term.label"])
    # print(len(clade_only))

    # Orig agg vs LP agg
    print("LP clade vs LP agg")
    both = pd.merge(lp_clade, lp_agg, how="inner", on="term.label")
    print(len(both))
    # print(both["term.label"])

    clade_only = pd.merge(orig_agg, both, how="left", on="term.label", indicator=True)
    clade_only = clade_only[clade_only["_merge"] == "left_only"]
    # print(orig_only["term.label"])
    # print(len(clade_only))

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
    vcs = all["term.label"].value_counts()
    print(vcs[vcs==3])

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
    compare_results(df)


if __name__ == "__main__":
    main()