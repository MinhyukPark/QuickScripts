import sys

from Bio import SeqIO
import dendropy

import decomposer

guide_tree = dendropy.Tree.get(path=sys.argv[1], schema="newick")
namespace = guide_tree.taxon_namespace
guide_tree.is_rooted = False
guide_tree.resolve_polytomies(limit=2)
guide_tree.collapse_basal_bifurcation()
guide_tree.update_bipartitions()

maxSubsetSize = int(sys.argv[2])
original_sequence_file = sys.argv[3]
output_partition_file_prefix = sys.argv[4]
trees = decomposer.decomposeTree(guide_tree, maxSubsetSize, mode="centroid")
clusters = []
for tree in trees:
    keep = [n.taxon.label for n in tree.leaf_nodes()]
    clusters.append(set(keep))
print(len(clusters))

files = [output_partition_file_prefix + str(i) + ".out" for i in range(len(clusters))]
sequence_partitions = [[] for _ in range(len(clusters))]

for sequence in SeqIO.parse(open(original_sequence_file), "fasta"):
    for cluster_index,cluster in enumerate(clusters):
        if(sequence.id in cluster):
            sequence_partitions[cluster_index].append(sequence)
for sequence_partition_index,sequence_partition in enumerate(sequence_partitions):
    SeqIO.write(sequence_partition, files[sequence_partition_index], "fasta")
