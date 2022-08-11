import sys

import click
import treeswift
from Bio import SeqIO


def resolve_polytomy_helper(input_tree, output_file):
    full_tree = treeswift.read_tree_newick(input_tree)
    full_tree.resolve_polytomies()
    if(hide_prefix):
        full_tree.write_tree_newick(output_file, hide_rooted_prefix=True)
    else:
        full_tree.write_tree_newick(output_file)

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--output-file", required=True, type=str, help="Output file path for the resolved tree")
def resolve_polytomy(input_tree, output_file):
    """This script resolve the polytomies in the input tree randomly with 0-length edges
    """
    resolve_polytomy_helper(input_tree, output_file)


if __name__ == "__main__":
    resolve_polytomy()
