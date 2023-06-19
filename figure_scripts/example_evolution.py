import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import nbinom
import argparse



data_params = {
        "model": "c",
        "n1": 15,
        "n2": 5,
        "t": 10,
        "theta_base": 100,
        "theta_ratio": 1.02,
        "alpha": 0.3,
        "sq": .25,
        "nb": False,
        "r": 100
}
plot_params = {
        "figsize":(7,4),
        "style1":{"color":"black", "linewidth":.6, "linestyle":"solid"},
        "style2":{"color":"black", "linewidth":.6, "linestyle":"solid"},
        # "ylim": {"y_min":142,"y_max":148,},
        "ylim": None,
        "title":"", 
        "save_path":None, 
        "save_type":"svg"
}

def get_bm(t, theta_base):
    w = np.random.default_rng().normal(0.0, scale=1.0, size=t)
    w = np.insert(w, 0, theta_base)
    return np.cumsum(w)

def get_ou(t, theta_base, alpha, theta, sq):
    w = np.random.default_rng().normal(0.0, scale=1.0, size=t)
    x = [theta_base]
    for i in range(t):
        dx = alpha * (theta - x[i]) + sq*w[i]
        x.append(x[i]+dx)
    return x

def add_nbinom(data, r):
    result = np.zeros((np.shape(data)[0], np.shape(data)[1]+1))
    for time_series in range(len(data)):
        p = r / (r + data[time_series][-1])
        result[time_series] = np.append(data[time_series], nbinom.rvs(n=r, p=p, size=1))
    return result

def get_save_name(data_params, plot_params):
    save_name = (plot_params["save_path"]
        + data_params["model"] +
        "_n1_" + str(data_params["n1"]) +
        "_n2_" + str(data_params["n2"]) + 
        "_t_" + str(data_params["t"]) +
        "_tb_" + str(data_params["theta_base"]) +
        "_tr_" + str(data_params["theta_ratio"]) +
        "_a_" + str(data_params["alpha"]) +
        "_sq_" + str(data_params["sq"])
    )
    if data_params["nb"]:
        save_name += "_r_" + str(data_params["r"])
    save_name += "." + plot_params["save_type"]
    return save_name

def make_plot(group1, group2, plot_params, save_name, show=False):
    # plt.figure()
    fig, axes = plt.subplots(nrows = 1, ncols = 2, sharey = True, 
                             gridspec_kw={'width_ratios': [.8,.2]},
                             figsize=plot_params["figsize"])
    for i in range(len(group1)):
        axes[0].plot(group1[i], c=plot_params["style1"]["color"], linewidth=plot_params["style1"]["linewidth"], linestyle=plot_params["style1"]["linestyle"])
    for i in range(len(group2)):
        axes[0].plot(group2[i], c=plot_params["style2"]["color"], linewidth=plot_params["style2"]["linewidth"], linestyle=plot_params["style2"]["linestyle"])
    axes[1].hist(group1[:,-1], orientation="horizontal", color=plot_params["style1"]["color"])
    axes[1].hist(group2[:,-1], orientation="horizontal", color=plot_params["style2"]["color"])
    if plot_params["ylim"]:
        axes[0].ylim(plot_params["ylim"]["y_min"], plot_params["ylim"]["y_max"])
    # fig.xticks([])
    # fig.yticks([])
    plt.title(plot_params["title"])
    # axes[1]
    # axes.xlabel("time")
    # axes.ylabel("expression value")
    if save_name:
        plt.savefig(save_name, format=plot_params["save_type"])
    if show:
        plt.show()

def get_data(model, n1, n2, t, theta_base, theta_ratio, alpha, sq, r, nb):
    theta_0 = (2*theta_base)/(1+theta_ratio)
    theta_1 = (2*theta_base*theta_ratio) / (1+theta_ratio)
    if model == "c":
        group1 = np.asarray([get_ou(t, theta_base, alpha, theta_1, sq) for i in range(n1)])
        group2 = np.asarray([get_ou(t, theta_base, alpha, theta_1, sq) for i in range(n2)])
    elif model == "a_a":
        group1 = np.asarray([get_ou(t, theta_base, alpha, theta_0, sq) for i in range(n1)])
        group2 = np.asarray([get_ou(t, theta_base, alpha, theta_1, sq) for i in range(n2)])
    elif model == "a_c":
        group1 = np.asarray([get_ou(t, theta_base, alpha, theta_base, sq) for i in range(n1)])
        group2 = np.asarray([get_ou(t, theta_base, alpha, theta_1, sq) for i in range(n2)])
    elif model == "a_n":
        group1 = np.asarray([get_bm(t, theta_base) for i in range(n1)])
        group2 = np.asarray([get_ou(t, theta_base, alpha, theta_0, sq) for i in range(n2)])
    else:
        group1 = np.asarray([get_bm(t, theta_base) for i in range(n1)])
        group2 = np.asarray([get_bm(t, theta_base) for i in range(n2)])
    if nb:
        group1 = add_nbinom(group1, r)
        group2 = add_nbinom(group2, r)
    return group1, group2
    
def main(data_params, plot_params, show=False):
    group1, group2 = get_data(
                        model = data_params["model"], 
                        n1 = data_params["n1"],
                        n2 = data_params["n2"],
                        t = data_params["t"],
                        theta_base = data_params["theta_base"],
                        theta_ratio = data_params["theta_ratio"],
                        alpha = data_params["alpha"],
                        sq = data_params["sq"],
                        r = data_params["r"],
                        nb = data_params["nb"])
    if plot_params["save_path"]:
        save_name = get_save_name(data_params, plot_params)
    else:
        save_name = None
    make_plot(group1, group2, plot_params, save_name, show)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--save_path", type=str, action="store", default=None)
    parser.add_argument("--title", type=str, action="store", default="")
    parser.add_argument("--show", "-s", action="store_true")
    args = parser.parse_args()
    plot_params["save_path"] = args.save_path
    plot_params["title"] = args.title
    main(data_params, plot_params, args.show)
