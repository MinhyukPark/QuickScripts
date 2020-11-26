import random
import sys

import click
import treeswift
from Bio import SeqIO

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--sequence-file", required=True, type=click.Path(exists=True), help="Sequence file with taxa labels to induce with")
@click.option("--output-prefix", required=True, type=str, help="Output file prefix for the induced subtree")
def split_testing_training(input_tree, sequence_file, output_prefix):
    """This script induces the input tree on the sequence file
    """
    full_tree = treeswift.read_tree_newick(input_tree)

    to_keep_node_labels = set()
    for sequence in SeqIO.parse(open(sequence_file), "fasta"):
        to_keep_node_labels.add(sequence.id)

    induced_tree = full_tree.extract_tree_with(to_keep_node_labels);

    induced_tree.write_tree_newick(output_prefix + "induced.tree")


if __name__ == "__main__":
    split_testing_training()
