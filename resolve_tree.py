import random
import sys

import click
import treeswift
import dendropy

from Bio import SeqIO

def resolve_tree_helper(input_tree, output_file, hide_prefix):
    # full_tree = treeswift.read_tree_newick(input_tree)
    # full_tree.resolve_polytomies()
    # full_tree.write_tree_newick(output_file)
    full_tree = dendropy.Tree.get(path=input_tree, schema="newick")
    full_tree.resolve_polytomies()
    # full_tree.is_rooted = False
    # full_tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    full_tree.update_bipartitions()
    full_tree.write(path=output_file, schema="newick", suppress_rooting=hide_prefix)

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--output-file", required=True, type=str, help="Output file path for the resolved tree")
@click.option("--hide-prefix", is_flag=True, help="Output file path for the resolved tree")
def resolve_tree(input_tree, output_file, hide_prefix):
    """This script resolves the input tree
    """
    resolve_tree_helper(input_tree, output_file, hide_prefix)


if __name__ == "__main__":
    resolve_tree()
