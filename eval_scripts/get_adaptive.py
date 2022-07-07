import pandas as pd
import os
import trisicell as tsc
import argparse


def main(result_path):
    sc_data = tsc.datasets.sublines_scrnaseq()
    exp_data = sc_data["expression"]

    results_df = pd.DataFrame()

    for f in os.listdir(result_path):
        temp_df = pd.read_csv(result_path+f)
        results_df = results_df.append(temp_df, ignore_index=True)

    results_df["ens_gene"] = exp_data.var.index[results_df["gene"].tolist()]

    print(results_df[results_df["adaptive"] == "adaptive"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("result_path", type=str, action="store")
    args = parser.parse_args()
    main(args.result_path)