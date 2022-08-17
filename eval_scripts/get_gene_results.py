from operator import mod
import pandas as pd
import os
import argparse


def main(result_path, save_path, model, inverse):
    if model == "c":
        model_col = "constrained_vs_neutral"
        cols = ["ensemble_id", "gene_name", "pvalue", "qvalue", model_col]
        if inverse:
            model_goal = "neutral"
        else:
            model_goal = "constrained"
    if model == "a":
        model_col = "adaptive"
        cols = ["ensemble_id", "gene_name","ou2_vs_bm_pvalue","ou2_vs_ou1_pvalue","qvalue","adaptive"]
        if inverse:
            model_goal = "not-adaptive"
        else:
            model_goal = "adaptive"

    results_df = pd.DataFrame()

    for f in [f for f in os.listdir(result_path) if ".csv" in f]:
        temp_df = pd.read_csv(result_path+f)
        results_df = results_df.append(temp_df, ignore_index=True)

    results_df[["ensemble_id", "gene_name"]] = results_df["gene"].str.split("_", expand=True)
    results_df = results_df.drop("gene", axis=1)
    results_df = results_df.reindex(columns=cols)
    pos_res = results_df[results_df[model_col] == model_goal]
    # print(pos_res)
    print(str(len(pos_res)) + " positive results out of " + str(len(results_df)))

    if save_path:
        # genes = pos_res["gene_name"].tolist()
        # with open(save_path, "w") as f:
        #     f.writelines("\n".join(genes))
        pos_res.to_csv(save_path)
        # temp = results_df[["ensemble_id", "gene_name"]]
        # temp.to_csv(save_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("result_path", type=str, action="store")
    parser.add_argument("model", type=str, action="store")
    parser.add_argument("-s", "--save_path", type=str, action="store", default=None)
    parser.add_argument("-i", "--inverse", action="store_true")
    args = parser.parse_args()
    if args.model not in ["c", "a"]:
        print("Model must be either 'c' for constrained or 'a' for adaptive.")
        exit()
    main(args.result_path, args.save_path, args.model, args.inverse)