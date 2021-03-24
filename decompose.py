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
@click.option("--fragmentary-sequence-file", required=False, type=click.Path(exists=True), help="Fragmentary aligned sequence file on some set of taxa")
@click.option("--output-prefix", required=True, type=str, help="Output file prefix for each subset")
@click.option("--maximum-size", required=True, type=int, help="Maximum size of output subsets")
@click.option("--longest-edge", is_flag=True, help="Specifying longest edge decomposition")
@click.option("--full-length", required=False, type=int, help="Specifyng a length for the sequences")
@click.option("--longest", required=False, is_flag=True, help="Samples from the longest sequences in the case of incomplete sampling")
@click.option("--incomplete", required=False, type=int, help="Specifying the size of incomplete decomposition (incomplete < maximum-size)")
@click.option("--support-threshold", required=False, type=float, help="Specifying a support threshold for decomposition", default=0.95)
@click.option("--mode", type=click.Choice(["centroid", "support", "heuristic"], case_sensitive=False), required=False, default="centroid")
def decompose_tree(input_tree, sequence_file, fragmentary_sequence_file, output_prefix, maximum_size, longest_edge, full_length, longest, incomplete, support_threshold, mode):
    '''This script decomposes the input tree and outputs induced alignments on the subsets.
    '''
    guide_tree = dendropy.Tree.get(path=input_tree, schema="newick")
    namespace = guide_tree.taxon_namespace
    guide_tree.is_rooted = False
    guide_tree.resolve_polytomies(limit=2)
    guide_tree.collapse_basal_bifurcation()
    guide_tree.update_bipartitions()

    # create mapping of closest sequencs
    fragmentary_mapping = {}
    for fragmentary_sequence in SeqIO.parse(open(fragmentary_sequence_file), "fasta"):
        current_id = fragmentary_sequence.id
        best_distance = len(fragmentary_sequence)
        best_id = None
        for full_length_sequence in SeqIO.parse(open(sequence_file), "fasta"):
            current_hamming_distance = 0
            for position,nucleotide in enumerate(fragmentary_sequence):
                if(fragmentary_sequence[position] != "-" and full_length_sequence[position] != "-"):
                    if(fragmentary_sequence[position] != full_length_sequence[position]):
                        current_hamming_distance += 1
            if(current_hamming_distance < best_distance):
                best_distance = current_hamming_distance
                best_id = full_length_sequence.id
        if(best_id not in fragmentary_mapping):
            fragmentary_mapping[best_id] = []
        fragmentary_mapping[best_id].append(current_id)

    trees = None
    sys.stderr.write(str(fragmentary_mapping) + "\n")
    if(fragmentary_sequence_file == None):
        trees = decomposer.decomposeTree(guide_tree, maximum_size, support_threshold, mode=mode)
    elif(fragmentary_sequence_file != None):
        trees = decomposer.fragmentary_decompose_tree(guide_tree, maximum_size, fragmentary_mapping)


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
                fragmentary_sequences = fragmentary_mapping[sequence.id]
                for fragmentary_sequence in fragmentary_sequences:
                    sequence_partitions[cluster_index].append(fragmentary_sequence)

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
