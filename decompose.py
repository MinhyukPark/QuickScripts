import random
import sys

import click
import dendropy
from Bio import SeqIO

import decomposer

# argv guide_tree, subset_size, alignment_on_full_taxa, output_file_prefix

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--sequence-file", required=True, type=click.Path(exists=True), help="Aligned sequence file on the full taxa")
@click.option("--output-prefix", required=True, type=str, help="Output file prefix for each subset")
@click.option("--maximum-size", required=True, type=int, help="Maximum size of output subsets")
@click.option("--longest-edge", is_flag=True, help="Specifying longest edge decomposition")
@click.option("--full-length", required=False, type=int, help="Specifyng a length for the sequences")
@click.option("--longest", required=False, is_flag=True, help="Samples from the longest sequences in the case of incomplete sampling")
@click.option("--incomplete", required=False, type=int, help="Specifying the size of incomplete decomposition (incomplete < maximum-size)")
@click.option("--support-threshold", required=False, type=float, help="Specifying a support threshold for decomposition", default=0.95)
def decompose_tree(input_tree, sequence_file, output_prefix, maximum_size, longest_edge, full_length, longest, incomplete, support_threshold):
    '''This script decomposes the input tree and outputs induced alignments on the subsets.
    '''
    guide_tree = dendropy.Tree.get(path=input_tree, schema="newick")
    namespace = guide_tree.taxon_namespace
    guide_tree.is_rooted = False
    guide_tree.resolve_polytomies(limit=2)
    guide_tree.collapse_basal_bifurcation()
    guide_tree.update_bipartitions()

    trees = decomposer.decomposeTree(guide_tree, maximum_size, support_threshold, mode="centroid")
    clusters = []
    for tree in trees:
        keep = [n.taxon.label.replace("_"," ") for n in tree.leaf_nodes()]
        clusters.append(set(keep))
    print(len(clusters))

    files = [output_prefix + str(i) + ".out" for i in range(len(clusters))]
    sequence_partitions = [[] for _ in range(len(clusters))]

    for sequence in SeqIO.parse(open(sequence_file), "fasta"):
        for cluster_index,cluster in enumerate(clusters):
            if(sequence.id.replace("_"," ") in cluster):
                sequence_partitions[cluster_index].append(sequence)

    for sequence_partition_index,sequence_partition in enumerate(sequence_partitions):
        SeqIO.write(sequence_partition, files[sequence_partition_index], "fasta")


    incomplete_sequences = []
    for sequence_partition_index,sequence_partition in enumerate(sequence_partitions):
        if(incomplete != None):
            if(longest):
                current_list = sorted(sequence_partitions[sequence_partition_index], key=lambda x:len(x.seq.ungap("-")), reverse=True)
                incomplete_sequences.extend(current_list[:incomplete])
            else:
                random.shuffle(sequence_partitions[sequence_partition_index])
                incomplete_sequences.extend(sequence_partitions[sequence_partition_index][:incomplete])
    if(incomplete != None):
        SeqIO.write(incomplete_sequences, output_prefix + "incomplete.out", "fasta")


if __name__ == "__main__":
    decompose_tree()
