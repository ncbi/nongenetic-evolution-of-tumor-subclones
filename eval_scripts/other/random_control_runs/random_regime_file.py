import csv
import dendropy
import random
from numpy.random import permutation
import argparse
from datetime import datetime
random.seed(datetime.now().timestamp())

def main(n):
    save_file = "regime_files/control/random{}.csv".format(n)
    tree = dendropy.Tree.get(path="tree_files/sc-bwes-cons-resolved-10.tree", schema="newick")
    chosen_nodes = []
    sample = permutation(tree.leaf_nodes())
    i=0
    while len(chosen_nodes) < k:
        if not sample[i].taxon.label in ["C1", "C4", "C22", "C11", "C15", "C16", "C18"]:
            chosen_nodes.append(sample[i])
        i+=1
    print([node.taxon.label for node in chosen_nodes])
    result = [["node", "node2", "regime"]]

    for leaf in tree.leaf_nodes():
        if leaf in chosen_nodes:
            result.append([leaf.taxon.label, "", "1_chosen"])
        else:
            result.append([leaf.taxon.label, "", "2_background"])

    internal_nodes = []
    for leaf1 in tree.leaf_nodes():
        for leaf2 in tree.leaf_nodes():
            if leaf1 == leaf2: 
                continue
            if leaf1.taxon.label == "C5" and leaf2.taxon.label == "C9": result.append(["C5", "C9", "2_background"]) # not sure why it skips this one
            anc = tree.mrca(taxa=[leaf1.taxon, leaf2.taxon])
            if anc in internal_nodes:
                continue
            result.append([leaf1.taxon.label, leaf2.taxon.label, "2_background"])
            internal_nodes.append(anc) 
    
    with open(save_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, action="store", help="n is the id of this run, only used in naming the file")
    args = parser.parse_args()
    main(args.n)