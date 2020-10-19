import random
import sys

import click
import treeswift
from Bio import SeqIO

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--sequence-file", required=True, type=click.Path(exists=True), help="Aligned sequence file on the full taxa")
@click.option("--output-prefix", required=True, type=str, help="Output file prefix for each subset")
@click.option("--minimum-size", required=False, type=int, help="Minimum size of each set")
def split_testing_training(input_tree, sequence_file, output_prefix, minimum_size):
    """This script decomposes the input tree and returns two clades by splitting at two random edges
    """
    guide_tree = treeswift.read_tree_newick(input_tree)
    num_nodes = guide_tree.num_nodes()
    testing_obtained = False
    training_obtained = False
    while(not testing_obtained or not training_obtained):
        testing_obtained = False
        training_obtained = False
        for node_index,node in enumerate(guide_tree.traverse_postorder(leaves=False, internal=True)):
            if(minimum_size < node_index):
                current_children = node.child_nodes()
                training_node = current_children[0]
                testing_node = current_children[1]
                training_tree = guide_tree.extract_subtree(training_node)
                testing_tree = guide_tree.extract_subtree(testing_node)
                if(training_tree.num_nodes() > minimum_size and testing_tree.num_nodes() > minimum_size):
                    testing_obtained = True
                    training_obtained = True
                    break

    testing_tree.write_tree_newick(output_prefix + "testing.tree")
    training_tree.write_tree_newick(output_prefix + "training.tree")

    current_cluster = [current_node.label.replace("_", " ") for current_node in testing_tree.traverse_postorder(leaves=True, internal=False)]
    current_sequence_partition = []
    for sequence in SeqIO.parse(open(sequence_file), "fasta"):
        if(sequence.id.replace("_", " ") in current_cluster):
            current_sequence_partition.append(sequence)
    SeqIO.write(current_sequence_partition, output_prefix + "testing.fasta", "fasta")

    current_cluster = [current_node.label.replace("_", " ") for current_node in training_tree.traverse_postorder(leaves=True, internal=False)]
    current_sequence_partition = []
    for sequence in SeqIO.parse(open(sequence_file), "fasta"):
        if(sequence.id.replace("_", " ") in current_cluster):
            current_sequence_partition.append(sequence)
    SeqIO.write(current_sequence_partition, output_prefix + "training.fasta", "fasta")

if __name__ == "__main__":
    split_testing_training()
