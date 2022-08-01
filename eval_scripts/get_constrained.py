import pandas as pd
import os
import trisicell as tsc
import argparse


def main(result_path, save_path):
    sc_data = tsc.datasets.sublines_scrnaseq()
    exp_data = sc_data["expression"]

    results_df = pd.DataFrame()

    for f in [f for f in os.listdir(result_path) if ".csv" in f]:
        temp_df = pd.read_csv(result_path+f)
        results_df = results_df.append(temp_df, ignore_index=True)

    results_df["ens_gene"] = exp_data.var.index[results_df["gene"].tolist()]

    constrained = results_df[results_df["constrained_vs_neutral"] == "constrained"]
    print(constrained)
    print(len(constrained))

    if save_path:
        genes = [g.split("_")[-1] for g in constrained["ens_gene"]]
        with open(save_path, "w") as f:
            f.writelines("\n".join(genes))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("result_path", type=str, action="store")
    parser.add_argument("-s", "--save_path", type=str, action="store", default=None)
    args = parser.parse_args()
    main(args.result_path, args.save_path)