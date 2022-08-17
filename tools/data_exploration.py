import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
import seaborn as sns
import os
import argparse

def violin_plot(df, save_path, show):
    fig = plt.figure(figsize=(10,5))
    sns.violinplot(x="species", y="exprval", data=df)
    if show:
        plt.show()
    if save_path:
        plt.savefig(save_path+"violin.png")

def get_subline_colors(df):
    green = ["C12", "C10", "C2", "C7", "C20"]
    purple = ["C17", "C3", "C9", "C24", "C8"]
    orange = ["C23", "C14", "C1", "C22", "C4"]
    blue = ["C5", "C6", "C19", "C21", "C13"]
    red = ["C15", "C18", "C11", "C16"]

    colors = []
    for subline in df["subline"]:
        if subline in green:
            colors.append("green")
        elif subline in purple:
            colors.append("purple")
        elif subline in orange:
            colors.append("orange")
        elif subline in blue:
            colors.append("blue")
        elif subline in red:
            colors.append("red")
        else: 
            colors.append("black")
    return colors

def get_agg_colors(df):
    colors = []
    for subline in df["subline"]:
        # if subline in ["C13", "C15", "C18", "C11", "C16", "C8", "C7", "C20"]:
        if subline in ["C18", "C15", "C11", "C16"]:
            colors.append("red")
        else:
            colors.append("black")
    return colors

def pca_plot(df, save_path, show):
    df["subclone"] = df["species"] + "_" + df["replicate"].astype(str)
    df = df.pivot_table(columns="gene", index="subclone", values="exprval")
    df = df.dropna(axis="index", how="all")
    df = df.fillna(0)
    
    pca_obj = PCA()
    pca_data = pca_obj.fit_transform(df)

    df = df.reset_index()
    df[["subline", "replicate"]] = df["subclone"].str.split("_", expand=True)
    df["pca1"] = pca_data[:, 0]
    df["pca2"] = pca_data[:, 1]
   
    df["subline_color"] = get_subline_colors(df)
    df["agg_color"] = get_agg_colors(df)
    
    fig = plt.figure()
    sns.scatterplot(x="pca1", y="pca2", data=df, hue="subline_color", palette=list(df["subline_color"].unique()))
    plt.title("PCA colored by subline")
    if show:
        plt.show()
    if save_path:
        plt.savefig(save_path+"pca_subline.png")

    fig = plt.figure()
    sns.scatterplot(x="pca1", y="pca2", data=df, hue="agg_color", palette=list(df["agg_color"].unique()))
    plt.title("PCA colored by most aggressive")
    if show:
        plt.show()
    if save_path:
        plt.savefig(save_path+"pca_agg.png")

def main(data_path, save_path, show):
    df = pd.DataFrame()
    for f in os.listdir(data_path):
        df = df.append(pd.read_csv(data_path+f))
    
    violin_plot(df, save_path, show)
    pca_plot(df, save_path, show)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data_path", type=str, action="store")
    parser.add_argument("--save_path", type=str, action="store", default=None)
    parser.add_argument("-v", "--show", action="store_true")
    args = parser.parse_args()
    if args.data_path[-1] != ["/"]: args.data_path += "/"
    main(args.data_path, args.save_path, args.show)