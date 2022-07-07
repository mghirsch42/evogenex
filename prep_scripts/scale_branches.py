from ast import Continue
import dendropy
from ete3 import PhyloTree
import numpy as np

tree_file = "tree_files/sc-bwes-cons.tree"
output_file = "tree_files/sc-bwes-cons-scaled-100.tree"

tree = dendropy.Tree.get(path=tree_file, schema="newick")
# t = PhyloTree(tree_file)
# t.show()

lineages = []
lineage_sums = {}
scaling_factors = {}

def rec(e, s, curr_lineage, lineages, lineage_sums):
    if e.length:
        s += e.length
    curr_lineage.append(e)
    lineages.append(curr_lineage)
    if e.is_leaf():
        lineage_sums[e] = s
        scaling_factors[e] = [e.length / s]
        return scaling_factors[e]
    for c in e.head_node.child_edges():
        sf = rec(c, s, curr_lineage, lineages, lineage_sums)
        if e in scaling_factors.keys():
            scaling_factors[e].extend([sf])
        else:
            scaling_factors[e] = [sf]
    return scaling_factors[e]

for e in tree.seed_node.child_edges():
    sf = rec(e, 0, [], lineages, lineage_sums)
    if e in scaling_factors.keys():
        scaling_factors[e].extend(sf)
    else:
        scaling_factors[e] = [sf]

# print(lineages)
# print(lineage_sums)
# print(scaling_factors)

def rec2(sfs):
    print("Enter:", sfs)
    print("ndim:", np.ndim(sfs), len(sfs))
    if not any(type(x) == list for x in sfs):
        print("Return:", sfs[0])
        return sfs[0]
    sfs2 = []
    for sub in sfs:
        sfs2.append(rec2(sub))
    print("Return:", np.mean(sfs2))
    return np.mean(sfs2)


for e in tree.edges():
    if not e in scaling_factors.keys(): continue
    print(scaling_factors[e])
    sf = rec2(scaling_factors[e])
    print(sf)
    e.length = e.length * sf

for e in tree.edges():
    e.length *= 100

print(tree.__str__())
t = PhyloTree(tree.__str__()+";")
t.show()

tree.write(path=output_file,
        schema="newick")