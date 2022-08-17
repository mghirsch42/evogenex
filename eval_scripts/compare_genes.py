from venv import create
import pandas as pd

def print_counts(orig_clade, orig_agg, lp_clade, lp_agg):
    print("Number of genes in each list")
    print("Orig, clade: {}".format(len(orig_clade)))
    print("Orig, agg: {}".format(len(orig_agg)))
    print("LP, clade: {}".format(len(lp_clade)))
    print("LP, agg: {}".format(len(lp_agg)))
    all = pd.concat([orig_clade, orig_agg, lp_clade, lp_agg])
    print("Total unique: {}".format(len(all["gene"].unique())))
    print()

def compare_trees(orig_clade, orig_agg, lp_clade, lp_agg):
    print("Orig clade vs LP clade")
    both = pd.merge(orig_clade, lp_clade, how="inner", on="gene")
    # print(both)
    print(len(both))

    print("Orig agg vs LP agg")
    both = pd.merge(orig_agg, lp_agg, how="inner", on="gene")
    # print(both)
    print(len(both))

def compare_adpt(orig_clade, orig_agg, lp_clade, lp_agg):
    print("Orig clade vs Orig agg")
    both = pd.merge(orig_clade, orig_agg, how="inner", on="gene")
    # print(both)
    print(len(both))

    print("LP clade vs LP agg")
    both = pd.merge(lp_clade, lp_agg, how="inner", on="gene")
    # print(both)
    print(len(both))

def compare_all(orig_clade, orig_agg, lp_clade, lp_agg):
    all = pd.concat([orig_clade, orig_agg, lp_clade, lp_agg])
    vcs = all["gene"].value_counts()
    print(vcs[vcs==4])

def create_table(orig_clade, orig_agg, lp_clade, lp_agg, save_path):
    all_res = pd.concat([orig_clade, orig_agg, lp_clade, lp_agg])
    all_genes = all_res[["ensemble_id", "gene_name"]]

    df = pd.DataFrame(columns=["ensemble_id", "gene_name", "orig_clade", "orig_agg", "scaled_clade", "scaled_agg"])
    df[["ensemble_id", "gene_name"]] = all_genes.drop_duplicates()
    df = df.reset_index(drop=True)
    df = df.fillna(0)

    orig_clade_genes = orig_clade["gene_name"].to_list()
    orig_agg_genes = orig_agg["gene_name"].to_list()
    lp_clade_genes = lp_clade["gene_name"].to_list()
    lp_agg_genes = lp_agg["gene_name"].to_list()

    df["orig_clade"] = [orig_clade[orig_clade["gene_name"] == g]["qvalue"].values[0] if g in orig_clade_genes else 0 for g in df["gene_name"] ]
    df["orig_agg"] = [orig_agg[orig_agg["gene_name"] == g]["qvalue"].values[0] if g in orig_agg_genes else 0 for g in df["gene_name"] ]
    df["lp_clade"] = [lp_clade[lp_clade["gene_name"] == g]["qvalue"].values[0] if g in lp_clade_genes else 0 for g in df["gene_name"] ]
    df["lp_agg"] = [lp_agg[lp_agg["gene_name"] == g]["qvalue"].values[0] if g in lp_agg_genes else 0 for g in df["gene_name"] ]

    df.to_csv(save_path, index=False)

def main():
    base_path = "results/tpm_allrep2/"
    orig_path = base_path + "orig/gene_lists/"
    lp_path = base_path + "lp_scaled/gene_lists/"

    orig_clade = pd.read_csv(orig_path+"adpt_clade/adpt_clade_gene_info.csv")
    # orig_clade.columns=["gene"]
    orig_agg = pd.read_csv(orig_path+"adpt_agg/adpt_agg_gene_info.csv")    
    # orig_agg.columns=["gene"]
    lp_clade = pd.read_csv(lp_path+"adpt_clade/adpt_clade_gene_info.csv")
    # lp_clade.columns=["gene"]
    lp_agg = pd.read_csv(lp_path+"adpt_agg/adpt_agg_gene_info.csv")
    # lp_agg.columns=["gene"]

    # print_counts(orig_clade, orig_agg, lp_clade, lp_agg)
    # compare_trees(orig_clade, orig_agg, lp_clade, lp_agg)
    # compare_adpt(orig_clade, orig_agg, lp_clade, lp_agg)
    # compare_all(orig_clade, orig_agg, lp_clade, lp_agg)
    create_table(orig_clade, orig_agg, lp_clade, lp_agg, "results/tpm_allrep2/gene_comp.csv")

if __name__ == "__main__":
    main()