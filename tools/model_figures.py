import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from math import ceil
from scipy.stats import nbinom
import argparse



def get_bm(t, theta_base, theta_ratio):

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
    result = data
    for time_series in range(len(data)):
        p = r / (r + data[time_series][-1])
        result[time_series][-1] = nbinom.rvs(n=r, p=p, size=1)
    return result

def plot_bm(bm1, bm2, label, save_name, show):
    plt.figure()
    for i in range(len(bm2)):
        plt.plot(bm2[i], c="black", linewidth=.6, linestyle="solid", label="Other sublines")
    for i in range(len(bm1)):
        plt.plot(bm1[i], c="red", linewidth=.6, linestyle=(0, (5, 5)), label=label)
    plt.ylim(140, 160)
    plt.title("Example Neutral Evolution")
    plt.xlabel("time")
    plt.ylabel("expression value")
    plt.savefig(save_name)
    if show:
        plt.show()

def plot_ou1(ou1a, ou1b, label, save_name, show):
    plt.figure()
    for i in range(len(ou1b)):
        plt.plot(ou1b[i], c="black", linewidth=.6, linestyle="solid", label="Other sublines")
    for i in range(len(ou1a)):
        plt.plot(ou1a[i], c="red", linewidth=.6, linestyle=(0, (5, 5)), label=label)
    plt.ylim(140, 160)
    plt.title("Example Constrained Evolution")
    plt.xlabel("time")
    plt.ylabel("expression value")
    plt.savefig(save_name)
    if show:
        plt.show()

def plot_ou2(ou1, ou2, label, save_name, show):
    plt.figure()
    for i in range(len(ou2)):
        plt.plot(ou2[i], c="black", linewidth=.6, linestyle="solid", label="Other sublines")
    for i in range(len(ou1)):
        plt.plot(ou1[i], c="red", linewidth=.6, linestyle=(0, (5, 5)), label=label)
    plt.ylim(140, 160)
    plt.title("Example Adaptive Evolution")
    plt.xlabel("time")
    plt.ylabel("expression value")
    plt.savefig(save_name)
    if show:
        plt.show()

def neutral():
    bm1 = np.asarray([get_bm(t, theta_base, theta_ratio) for i in range(9)])
    bm2 = np.asarray([get_bm(t, theta_base, theta_ratio) for i in range(14)])
    print(bm1)
    if nb == "y":
        bm1 = add_nbinom(bm1, r)
        bm2 = add_nbinom(bm2, r)
    plot_bm(bm1, bm2, label="Aggressive Responder 8", 
                    save_name="figures/models/neut_an8.png", 
                    show=True)

# Constrained Evolution
def constrained():
    ou1a = np.asarray([get_ou(t, theta_base, alpha, theta_1, sq) for i in range(9)])
    ou1b = np.asarray([get_ou(t, theta_base, alpha, theta_1, sq) for i in range(14)])
    if nb == "y":
        ou1a = add_nbinom(ou1a, r)
        ou1b = add_nbinom(ou1b, r)
    plot_ou1(ou1a, ou1b, label="Aggressive Responder 8", 
                    save_name="figures/models/const_an8_a_{}_sq_{}_tr_{}_r_{}_nb_.png".format(alpha, sq, theta_base, r, nb), 
                    show=True)

# Adaptive Evolution
def adaptive():
    ou1 = np.asarray([get_ou(t, theta_base, alpha, theta_0, sq) for i in range(9)])
    ou2 = np.asarray([get_ou(t, theta_base, alpha, theta_1, sq) for i in range(14)])
    if nb == "y":
        ou1 = add_nbinom(ou1, r)
        ou2 = add_nbinom(ou2, r)
    plot_ou2(ou1, ou2, label="Aggressive Responder 8", 
                    save_name="figures/models/adpt_an8_a_{}_sq_{}_tr_{}_r_{}_nb_.png".format(alpha, sq, theta_base, r, nb), 
                    show=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("models", type=str, action="store", nargs="+")
    args = parser.parse_args()

    plt.rcParams["figure.figsize"] = (6,4)
    t = 10
    n = 23
    alpha = 1
    sq = 1.5
    theta_base = 150
    theta_ratio = 1.04
    theta_0 = (2*theta_base)/(1+theta_ratio)
    theta_1 = (2*theta_base*theta_ratio) / (1+theta_ratio)
    r=1000
    nb = "n"

    if "n" in args.models:
        neutral()
    if "c" in args.models:
        constrained()
    if "a" in args.models:
        adaptive()

