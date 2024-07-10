import pandas as pd
import os
# from matplotlib import pyplot as plt
# import seaborn as sns

def get_neut_sim_data(results_folder, regime, col):
    df = pd.DataFrame(columns=["sigma_sq", "r"])
    for f in os.listdir(results_folder):
        temp = pd.read_csv(results_folder+f)
        if regime == "adpt":
            temp = temp.iloc[::2, :] # Deal with double rows with the different theta values
            temp = temp.reset_index()
        temp["sigma_sq"] = f.split("_")[1]
        temp["r"] = f.split("_")[3][:-4]
        temp = temp[["sigma_sq", "r", col]]
        df = pd.concat([df, temp], ignore_index=True)
    return df

def get_const_sim_data(results_folder, regime, col):
    df = pd.DataFrame(columns=["sigma_sq", "r", "alpha"])
    for f in os.listdir(results_folder):
        temp = pd.read_csv(results_folder+f)
        if regime == "adpt":
            temp = temp.iloc[::2, :] # Deal with double rows with the different theta values
            temp = temp.reset_index()
        temp["sigma_sq"] = f.split("_")[3]
        temp["r"] = f.split("_")[5]
        temp["alpha"] = f.split("_")[7][:-4]
        temp = temp[["sigma_sq", "r", "alpha", col]]
        df = pd.concat([df, temp], ignore_index=True)
    return df

def get_adpt_sim_data(results_folder, regime, col):
    df = pd.DataFrame(columns=["sigma_sq", "r", "alpha", "theta_ratio"])
    for f in os.listdir(results_folder):
        temp = pd.read_csv(results_folder+f)
        if regime == "adpt":
            temp = temp.iloc[::2, :] # Deal with double rows with the different theta values
            temp = temp.reset_index()
        temp["theta_ratio"] = f.split("_")[1]
        temp["sigma_sq"] = f.split("_")[3]
        temp["r"] = f.split("_")[5]
        temp["alpha"] = f.split("_")[7][:-4]
        temp = temp[["sigma_sq", "r", "alpha", "theta_ratio", col]]
        df = pd.concat([df, temp], ignore_index=True)
    return df

def group(df, sim, regime, sim_vars, result_col, result_goal):
    df["{} sim, {} pred".format(sim, regime)] = (df[result_col] == result_goal).astype(int)
    count = df.groupby(sim_vars).count()
    gb = df.groupby(sim_vars).sum().reset_index()
    return gb

def main(sims, regimes, results_folder, save_folder, sim_var_dict):
    for sim_data in sims:
        print(sim_data)
        for regime in regimes:
            print("--", regime)
            if regime == "adpt":
                col = "adaptive"
                goal = "adaptive"
            else:
                col = "constrained_vs_neutral"
                goal = "constrained"
            if sim_data == "neut":
                df = get_neut_sim_data(results_folder.format(sim_data, regime), regime, col)
                gb = group(df, sim_data, regime, sim_var_dict[sim_data], col, goal)
            elif sim_data == "adpt":
                df = get_adpt_sim_data(results_folder.format(sim_data, regime), regime, col)
                gb = group(df, sim_data, regime, sim_var_dict[sim_data], col, goal)
            elif sim_data == "const":
                df = get_const_sim_data(results_folder.format(sim_data, regime), regime, col)
                gb = group(df, sim_data, regime, sim_var_dict[sim_data], col, goal)
            else:
                print("sim_data not recognized, exiting") 
                exit()
            gb.to_csv(save_folder.format(sim_data, regime), index=False)

if __name__ == "__main__":
    sims = ["neut", "const", "adpt"]
    regimes = ["const", "adpt"]
    results_folder = "results/simulated/{}_sim_{}_run/"
    save_folder = "results/simulated/summaries/{}_sim_{}_run.csv"

    sim_var_dict = {
        "neut": ["sigma_sq", "r"],
        "const": ["sigma_sq", "r", "alpha"],
        "adpt": ["sigma_sq", "r", "alpha", "theta_ratio"]
    }

    main(sims, regimes, results_folder, save_folder, sim_var_dict)