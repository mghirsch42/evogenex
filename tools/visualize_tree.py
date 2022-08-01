from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node
import pandas as pd
import argparse


def main(newick_file, regime_file, show, output_file):

    t = Tree(newick_file)

    if regime_file:
        regimes = pd.read_csv(regime_file)
        for idx, row in regimes.iterrows():
            if pd.isna(row["node2"]):
                n = t.get_leaves_by_name(row["node"])[0]
                f = TextFace(row["node"])
                f.rotation = -90
                n.add_face(f, column=0)
            else:
                n = t.get_common_ancestor(row["node"], row["node2"])
            mode = row["regime"]
            ns = NodeStyle()
            ns["size"] = 12
            if mode == "aggressive":
                ns["fgcolor"] = "red"
            elif mode == "partaggressive":
                ns["fgcolor"] = "blue"
            else:
                ns["fgcolor"] = "black"
            n.set_style(ns)
    
    t.get_children()[0].dist = 0  # Set the root branch to have 0 length

    ts = TreeStyle()
    ts.show_branch_length = True
    ts.rotation = 90
    ts.show_scale = False
    ts.show_leaf_name = False

    if output_file:
        t.render(output_file, tree_style=ts)

    if show:
        t.show(tree_style=ts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("newick_file", type=str, action="store")
    parser.add_argument("--regime_file", type=str, action="store", default=None)
    parser.add_argument("-s", "--show", action="store_true")
    parser.add_argument("--output_file", type=str, action="store", default=None)
    args = parser.parse_args()
    main(args.newick_file, args.regime_file, args.show, args.output_file)