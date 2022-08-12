from ast import arg
from venv import create
import pandas as pd
import argparse

def create_subtable(go_file, ann_data_set):
    # with open(go_file, "r") as f:
    #     lines = f.readlines()
    # data = list(map(str.split, map(str.strip, lines), ["\t"]*len(lines)))
    # raw_df = pd.DataFrame(data[1:], columns=["Description (GO)", "Ref list", "upload", "expected", "over/under", "fold enrichment", "pvalue", "FDR"])
    raw_df = pd.read_csv(go_file)
    # print(raw_df)
    latex_df = pd.DataFrame(columns=["GO ID", "Description", "Fold enrichment", "P-value", "FDR"])
    latex_df["GO ID"] = raw_df["term.id"]
    latex_df["Description"] = raw_df["term.label"]
    latex_df["Fold enrichment"] = "$" + raw_df["plus_minus"] + raw_df["fold_enrichment"].round(2).astype(str) + "$"
    pvals = ["{:.2e}".format(p).split("e") for p in raw_df["pValue"]]
    pvals = [[float(p[0]), int(p[1])] for p in pvals]
    latex_df["P-value"] = ["${}\\times10^{{{}}}$".format(p[0], p[1]) for p in pvals]
    
    fdr = ["{:.2e}".format(f).split("e") for f in raw_df["fdr"]]
    fdr = [[float(f[0]), int(f[1])] for f in fdr]
    latex_df["FDR"] = ["${}\\times10^{{{}}}$".format(f[0], f[1]) for f in fdr]
    
    table_str = ""
    if len(latex_df) > 1:
        table_str += "\\multirow[t]{" + str(len(latex_df)) + "}{\\linewidth}{" + "\\\\".join(ann_data_set.split(" ")) + "}\n"
    else:
        table_str += "\\\\".join(ann_data_set.split(" ")) + "\n"
    for idx, row in latex_df.iterrows():
        table_str += "& " + " & ".join(row) + " \\\\ \\cline{2-6}\n"
    table_str = table_str[:-12] # remove last cline
    table_str += "\\hline\n"
    return table_str

def main(input_files, ann_data_sets, save_path, verbose):
    table_str = ("\\begin{longtblr}[\n" +
                 "caption={},\n" +
                 "label={}\n" + 
                 "]{|p{.14\\textwidth}|p{.16\\textwidth}|p{.3\\textwidth}|p{.1\\textwidth}|p{.15\\textwidth}|p{.15\\textwidth}|}\n"
                 "\\hline\n" +
                 "Annotation Data Set & ID & Description & Fold Enrich-ment & P-value & FDR \\\\ \\hline\n")

    for i in range(len(input_files)):
        table_str += create_subtable(input_files[i], ann_data_sets[i])

    table_str += ("\\end{longtblr}\n")
    
    if verbose:
        print(table_str)

    if save_path:
        with open(save_path, "w") as f:
            f.writelines(table_str)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_files", type=str, nargs="+", action="store", required=True)
    parser.add_argument("-a", "--ann_data_sets", type=str, nargs="+", action="store", required=True)
    parser.add_argument("-s", "--save_path", type=str, action="store", default=None)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    if len(args.input_files) != len(args.ann_data_sets):
        print("There must be the same numer of input files as ann data set labels.")
        exit()
    main(args.input_files, args.ann_data_sets, args.save_path, args.verbose)