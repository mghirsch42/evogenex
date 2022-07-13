import dendropy
from ete3 import PhyloTree
from matplotlib.pyplot import sca
import numpy as np
import argparse

    
def get_lineage_lengths_rec(edge, curr_sum, curr_lineage, lineage_lengths):
    if edge.length:
        curr_sum += edge.length
    curr_lineage.append(edge)
    if edge.is_leaf():
        lineage_lengths[edge] = curr_sum
    else:
        for child in edge.head_node.child_edges():
            get_lineage_lengths_rec(child, curr_sum, curr_lineage.copy(), lineage_lengths)

def get_scaling_factors_rec(edge, scaling_factors):
    if edge.is_leaf():
        return scaling_factors[edge]
    child_sfs = []
    for child_edge in edge.head_node.child_edges():
        child_sfs.append(get_scaling_factors_rec(child_edge, scaling_factors))
    scaling_factors[edge] = np.mean(child_sfs)
    return scaling_factors[edge]

def main(tree_file, output_file):
    
    # Read in tree
    tree = dendropy.Tree.get(path=tree_file, schema="newick")

    # Initialize dictionaries
    lineage_lengths = {} # Dictionary with terminal edges as keys and root to tip branch lengths for that edge as values
    scaling_factors = {} # Dictionary with edges as keys and float scaling factors as values

    # Calculate the lengths of the lineages and populate lineage_lengths
    for edge in tree.seed_node.child_edges():
        get_lineage_lengths_rec(edge, 0, [], lineage_lengths)
    min_len = min(lineage_lengths.values())
    print(min_len)
    # Calculate scaling factors of the terminal edges and populate those in scaling_factors
    for terminal_edge in lineage_lengths.keys():
        sf = min_len / lineage_lengths[terminal_edge]
        scaling_factors[terminal_edge] = sf

    # Calculate scaling factors for the rest of the edges and populate scaling_factors
    for edge in tree.seed_node.child_edges():
        get_scaling_factors_rec(edge, scaling_factors)

    # Multiply edges by scaling factors
    for edge in tree.edges():
        if edge in scaling_factors.keys():  # This will include a root edge that isn't scaled
            print("old length:", edge.length)
            print("scaling factor:", scaling_factors[edge])
            edge.length *= scaling_factors[edge]
            print("new length", edge.length)

    # Write to output
    tree.write(path=output_file, schema="newick")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree_file", type=str, action="store", default="tree_files/sc-bwes-cons.tree")
    parser.add_argument("--output_file", type=str, action="store", default="tree_files/sc-bwes-cons-scaled.tree")
    args = parser.parse_args()
    main(args.tree_file, args.output_file)