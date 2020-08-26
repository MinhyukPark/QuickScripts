import sys

import dendropy

def preprocess_tree(tree):
    tree.is_rooted = False
    tree.resolve_polytomies(limit=2)
    tree.collapse_basal_bifurcation()
    tree.update_bipartitions(suppress_unifurcations=False)

    return tree

left_tree_path = sys.argv[1]
right_tree_path = sys.argv[2]
left_tree_test = dendropy.Tree.get(path=left_tree_path, schema="newick")
right_tree_test = dendropy.Tree.get(path=right_tree_path, schema="newick")
left_tree_taxa_list = sorted([str(taxon) for taxon in left_tree_test.taxon_namespace])
right_tree_taxa_list = sorted([str(taxon) for taxon in right_tree_test.taxon_namespace])
if(left_tree_taxa_list != right_tree_taxa_list):
    print("ERROR: tree taxa lists differ")
    exit(1)


dataset = dendropy.DataSet()
namespace = dataset.new_taxon_namespace()
dataset.read(path=left_tree_path, schema="newick", taxon_namespace=namespace)
dataset.read(path=right_tree_path, schema="newick", taxon_namespace=namespace)
left_tree = dataset.tree_lists[0][0]
right_tree = dataset.tree_lists[1][0]

left_tree = preprocess_tree(left_tree)
right_tree = preprocess_tree(right_tree)

print(f"Num Taxa: {len(namespace)}")
unweighted_rf_distance = dendropy.calculate.treecompare.unweighted_robinson_foulds_distance(left_tree, right_tree)
rf_rate = unweighted_rf_distance / ((2 * len(namespace)) - 6)
print(f"unweighted_rf_distance: {unweighted_rf_distance}")
print(f"rf rate: {rf_rate}")
