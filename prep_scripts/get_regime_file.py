import csv
import dendropy

result = [["node", "node2", "regime"]]

tree = dendropy.Tree.get(path="sc-bwes-cons.tree", schema="newick")

for leaf in tree.leaf_nodes():
    result.append([leaf.taxon.label, "", "global"])

internal_nodes = []
for leaf1 in tree.leaf_nodes():
    for leaf2 in tree.leaf_nodes():
        if leaf1 == leaf2: 
            continue
        anc = tree.mrca(taxa=[leaf1.taxon, leaf2.taxon])
        if anc in internal_nodes: 
            continue
        result.append([leaf1.taxon.label, leaf2.taxon.label, "global"])
        internal_nodes.append(anc)
        
with open("test_regime.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(result)