import random
import sys

import click
import treeswift
from Bio import SeqIO


def induce_tree_helper(input_tree, sequence_file, output_file, hide_prefix):
    full_tree = treeswift.read_tree_newick(input_tree)

    to_keep_node_labels = set()
    for sequence in SeqIO.parse(open(sequence_file), "fasta"):
        to_keep_node_labels.add(sequence.id)

    induced_tree = full_tree.extract_tree_with(to_keep_node_labels);
    if(hide_prefix):
        induced_tree.write_tree_newick(output_file, hide_rooted_prefix=True)
    else:
        induced_tree.write_tree_newick(output_file)

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--sequence-file", required=True, type=click.Path(exists=True), help="Sequence file with taxa labels to induce with")
@click.option("--output-file", required=True, type=str, help="Output file path for the induced subtree")
@click.option("--hide-prefix", required=False, is_flag=True, help="Whether to include the rooted prefix in the tree file or not")
def induce_tree(input_tree, sequence_file, output_file, hide_prefix):
    """This script induces the input tree on the sequence file
    """
    induce_tree_helper(input_tree, sequence_file, output_file, hide_prefix)


if __name__ == "__main__":
    induce_tree()
