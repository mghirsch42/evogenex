import csv
import dendropy

result = [["node", "node2", "regime"]]

tree = dendropy.Tree.get(path="tree_files/resolved/clade_sub/sc-bwes-cons-resolved-10-clade.tree", schema="newick")

for leaf in tree.leaf_nodes():
    result.append([leaf.taxon.label, "", "global"])

internal_nodes = []
for leaf1 in tree.leaf_nodes():
    for leaf2 in tree.leaf_nodes():
        print(leaf1.taxon.label, leaf2.taxon.label)
        # if leaf1.taxon.label == "C5" or leaf2.taxon.label == "C5": print(leaf1.taxon.label, leaf2.taxon.label)
        if leaf1 == leaf2: 
            continue
        anc = tree.mrca(taxa=[leaf1.taxon, leaf2.taxon])
        print(anc)
        if anc in internal_nodes:
            continue
        print("using")
        result.append([leaf1.taxon.label, leaf2.taxon.label, "global"])
        internal_nodes.append(anc)

with open("regime_files/clade_sub/single_resolved.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(result)