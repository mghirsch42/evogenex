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
        cols = ["ensemble_id", "gene_name"]
        if inverse:
            model_goal = "not-adaptive"
        else:
            model_goal = "adaptive"

    results_df = pd.DataFrame(columns=cols)

    for f in [f for f in os.listdir(result_path) if ".csv" in f]:
        temp_df = pd.read_csv(result_path+f)
        results_df = results_df.append(temp_df, ignore_index=True)

    results_df[["ensemble_id", "gene_name"]] = results_df["gene"].str.split("_", expand=True)
    results_df = results_df.drop("gene", axis=1)
    
    df2 = pd.DataFrame(columns=["ensemble_id", "gene_name","ou1_conv","ou1_theta","ou1_alpha","ou1_sigma_sq","ou1_gamma","ou1_loglik","ou2_conv","ou2_theta","ou2_theta_base","theta_diff","ou2_alpha","ou2_sigma_sq","ou2_gamma","ou2_loglik","brown_conv","brown_theta","brown_sigma_sq","brown_gamma","brown_loglik","ou2_vs_bm_pvalue","ou2_vs_ou1_pvalue","adaptive"])

    for g in results_df["gene_name"].unique():
        sub = results_df[results_df["gene_name"] == g]  # Two entries - adaptive, then base
        df2 = df2.append(sub.iloc[0]) # Add adaptive entry
        df2.loc[df2["gene_name"] == g, "ou2_theta_base"] = sub["ou2_theta"].iloc[1] # Add base value

    df2["theta_diff"] = df2["ou2_theta"] - df2["ou2_theta_base"]

    # results_df = results_df.reindex(columns=cols)
    pos_res = df2[df2[model_col] == model_goal]
    pos_res = pos_res.drop(model_col, axis=1)
    # print(pos_res)
    print(str(len(pos_res)) + " positive results out of " + str(len(df2)))

    if save_path:
        pos_res.to_csv(save_path + "gene_info.csv", index=False)
        # with open(save_path+"genes.csv", "w") as f:
        #     f.writelines("\n".join(pos_res["gene_name"].unique()))
        
        # results_df.to_csv(save_path + "gene_info.csv", index=False)
        # with open(save_path+"genes.csv", "w") as f:
        #     f.writelines("\n".join(results_df["gene_name"].unique()))

    else:
        print(pos_res)

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