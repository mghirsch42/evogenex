import pandas as pd
import os
import argparse


def main(result_path, save_path, model):
    if model == "c":
        model_col = "constrained_vs_neutral"
        model_goal = "constrained"
    if model == "a":
        model_col = "adaptive"
        model_goal = "adaptive"

    results_df = pd.DataFrame()

    for f in [f for f in os.listdir(result_path) if ".csv" in f]:
        temp_df = pd.read_csv(result_path+f)
        results_df = results_df.append(temp_df, ignore_index=True)

    results_df[["ensemble_id", "gene_name"]] = results_df["gene"].str.split("_", expand=True)
    results_df = results_df.drop("gene", axis=1)
    results_df = results_df.reindex(columns=["ensemble_id", "gene_name", "pvalue", "qvalue", model_col])
    pos_res = results_df[results_df[model_col] == model_goal]
    # print(pos_res)
    print(str(len(pos_res)) + " positive results out of " + str(len(results_df)))

    if save_path:
        genes = pos_res["gene_name"].tolist()
        with open(save_path, "w") as f:
            f.writelines("\n".join(genes))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("result_path", type=str, action="store")
    parser.add_argument("model", type=str, action="store")
    parser.add_argument("-s", "--save_path", type=str, action="store", default=None)
    args = parser.parse_args()
    if args.model not in ["c", "a"]:
        print("Model must be either 'c' for constrained or 'a' for adaptive.")
        exit()
    main(args.result_path, args.save_path, args.model)