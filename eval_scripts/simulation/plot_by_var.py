import pandas as pd
import os
from matplotlib import pyplot as plt
import seaborn as sns

sim_data = "neut" # neut, const, or agg
regime = "const" # const or agg
results_folder = "results/simulated3/{}_sim_{}_run/".format(sim_data, regime)
save_folder = "figures/simulation3/{}_sim_{}_run/".format(sim_data, regime)
# groupby_var = "sigma_sq"
# groupby_var = "alpha"
groupby_var = "r"
# groupby_var = "theta_ratio"

def plot_by_var(df, var, result_type, title, save_file=None):
    sns.boxplot(data=df, x=var, y=result_type)
    plt.title(title)
    if save_file:
        plt.savefig(save_file)
    plt.show()

if sim_data == "neut":
    df = pd.DataFrame(columns=["sigma_sq", "r"])
    for f in os.listdir(results_folder):
        temp = pd.read_csv(results_folder+f)
        if regime == "agg":
            temp = temp.iloc[::2, :] # Deal with double rows with the different theta values
            temp = temp.reset_index()
        temp["sigma_sq"] = f.split("_")[1]
        temp["r"] = f.split("_")[3][:-4]
        df = pd.concat([df, temp], ignore_index=True)
    if regime == "const":
        df["neut sim, const pred"] = (df["constrained_vs_neutral"] == "const").astype(int)
        gb = df.groupby(["sigma_sq", "r"]).sum().reset_index()
        plot_by_var(gb, groupby_var, "neut sim, const pred", groupby_var+" neut sim, const pred\n(FP-constrained on neutral data)", save_file=save_folder+groupby_var+".png")
    if regime == "agg":
        df["neut sim, adpt pred"] = (df["adaptive"] == "adaptive").astype(int)
        gb = df.groupby(["sigma_sq", "r"]).sum().reset_index()
        plot_by_var(gb, groupby_var, "neut sim, adpt pred", groupby_var+" neut sim, adpt prep\n(FP-adaptive on neutral data)", save_file=save_folder+groupby_var+".png")
        print(df["neut sim, adpt pred"].sum())
    print(df)
    print(gb)
elif sim_data == "agg":
    df = pd.DataFrame(columns=["theta_ratio", "sigma_sq", "r", "alpha"])
    for f in os.listdir(results_folder):
        temp = pd.read_csv(results_folder+f)
        if regime == "agg":
            temp = temp.iloc[::2, :] # Deal with double rows with the different theta values
            temp = temp.reset_index()
        temp["theta_ratio"] = f.split("_")[1]
        temp["sigma_sq"] = f.split("_")[3]
        temp["r"] = f.split("_")[5]
        temp["alpha"] = f.split("_")[7][:-4]
        df = pd.concat([df, temp], ignore_index=True)
    print(df)
    if regime == "agg":
        df["adpt sim, adpt pred"] = (df["adaptive"] == "adaptive").astype(int)
        gb = df.groupby(["theta_ratio", "sigma_sq", "r", "alpha"]).sum().reset_index()
        plot_by_var(gb, groupby_var, "adpt sim, adpt pred", groupby_var+" adpt sim, adpt pred\n(TP-adaptive on adaptive data)", save_file=save_folder+groupby_var+".png")
    else:
        df["adpt sim, const pred"] = (df["constrained_vs_neutral"] == "constrained").astype(int)
        gb = df.groupby(["theta_ratio", "sigma_sq", "r", "alpha"]).sum().reset_index()
        plot_by_var(gb, groupby_var, "adpt sim, const pred", groupby_var+" adpt sim, const pred\n(FP-constrained on adaptive data)", save_file=save_folder+groupby_var+".png")
    print(gb)
elif sim_data == "const":
    df = pd.DataFrame(columns=["sigma_sq", "r", "alpha"])
    for f in os.listdir(results_folder):
        temp = pd.read_csv(results_folder+f)
        if regime == "agg":
            temp = temp.iloc[::2, :] # Deal with double rows with the different theta values
            temp = temp.reset_index()
        temp["sigma_sq"] = f.split("_")[3]
        temp["r"] = f.split("_")[5]
        temp["alpha"] = f.split("_")[7][:-4]
        df = pd.concat([df, temp], ignore_index=True)
    if regime == "const":
        df["const sim, const pred"] = (df["constrained_vs_neutral"] == "const").astype(int)
        gb = df.groupby(["sigma_sq", "r", "alpha"]).sum().reset_index()
        plot_by_var(gb, groupby_var, "const sim, const pred", groupby_var+" const sim, const pred\n(TP-constrained on constrained data)", save_file=save_folder+groupby_var+".png")
    if regime == "agg":
        df["const sim, adpt pred"] = (df["adaptive"] == "adaptive").astype(int)
        gb = df.groupby(["sigma_sq", "r", "alpha"]).sum().reset_index()
        plot_by_var(gb, groupby_var, "const sim, adpt pred", groupby_var+" const sim, adpt pred\n(FP-adaptive on constrained data)", save_file=save_folder+groupby_var+".png")
    print(df)
    print(gb.sum())
else:
    print("sim_data not recognized, exiting") 
    exit()

# print(df)