import random
import sys

import click
import treeswift
from Bio import SeqIO

def induce_tree_helper(input_tree, input_data, output_file, hide_prefix, input_type, resolve_polytomies):
    full_tree = treeswift.read_tree_newick(input_tree)

    to_keep_node_labels = set()
    if input_type == "fasta":
        for sequence in SeqIO.parse(open(input_data), "fasta"):
            to_keep_node_labels.add(sequence.id)
    elif input_type == "newick":
        to_keep_tree = treeswift.read_tree_newick(input_data)
        for current_node in to_keep_tree.traverse_leaves():
            to_keep_node_labels.add(current_node.label)


    induced_tree = full_tree.extract_tree_with(to_keep_node_labels);
    if(resolve_polytomies):
        induced_tree.resolve_polytomies()
    if(hide_prefix):
        induced_tree.write_tree_newick(output_file, hide_rooted_prefix=True)
    else:
        induced_tree.write_tree_newick(output_file)

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--sequence-file", required=True, type=click.Path(exists=True), help="Sequence file with taxa labels to induce with")
@click.option("--output-file", required=True, type=str, help="Output file path for the induced subtree")
@click.option("--input-type", type=click.Choice(["fasta", "newick"], case_sensitive=False), required=False, help="set it to newick if the input sequence file is actually a tree", default="fasta")
@click.option("--hide-prefix", required=False, is_flag=True, help="Whether to include the rooted prefix in the tree file or not", default=True)
@click.option("--resolve-polytomies", required=False, is_flag=True, help="Whether to randomly resolve polytomies with 0-length edges")
def induce_tree(input_tree, sequence_file, output_file, input_type, hide_prefix, resolve_polytomies):
    """This script induces the input tree on the sequence file
    """
    induce_tree_helper(input_tree, sequence_file, output_file, hide_prefix, input_type, resolve_polytomies)


if __name__ == "__main__":
    induce_tree()
