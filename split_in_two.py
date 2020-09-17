import sys

import click
import dendropy
from Bio import SeqIO

import decomposer

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--sequence-file", required=True, type=click.Path(exists=True), help="Aligned sequence file on the full taxa")
@click.option("--output-prefix", required=True, type=str, help="Output file prefix for each subset")
@click.option("--minimum-size", required=False, type=int, help="Minimum size of input alignment")
def split_in_two(input_tree, sequence_file, output_prefix, minimum_size):
    """This script decomposes the input tree into two and outputs induced alignments on the subsets
    """
    guide_tree = dendropy.Tree.get(path=input_tree, schema="newick")
    namespace = guide_tree.taxon_namespace
    guide_tree.is_rooted = False
    guide_tree.resolve_polytomies(limit=2)
    guide_tree.collapse_basal_bifurcation()
    guide_tree.update_bipartitions()

    centroid_edge = decomposer.getCentroidEdge(guide_tree)
    tree1,tree2 = decomposer.bipartitionByEdge(guide_tree, centroid_edge)
    trees = [tree1, tree2]

    if(len(guide_tree.leaf_nodes()) < minimum_size):
        print("-1")
        return -1

    clusters = []
    for tree in trees:
        keep = [n.taxon.label for n in tree.leaf_nodes()]
        clusters.append(set(keep))
    print(len(clusters))

    files = [output_prefix + "L.out", output_prefix + "R.out"]
    sequence_partitions = [[] for _ in range(len(clusters))]

    for sequence in SeqIO.parse(open(sequence_file), "fasta"):
        for cluster_index,cluster in enumerate(clusters):
            if(sequence.id in cluster):
                sequence_partitions[cluster_index].append(sequence)
    for sequence_partition_index,sequence_partition in enumerate(sequence_partitions):
        SeqIO.write(sequence_partition, files[sequence_partition_index], "fasta")

if __name__ == "__main__":
    split_in_two()
