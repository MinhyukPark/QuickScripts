import sys

import click
import treeswift
from Bio import SeqIO


def relabel_tree_helper(input_tree, tax_list, output_file, hide_prefix):
    full_tree = treeswift.read_tree_newick(input_tree)
    tax_map = {}
    with open(tax_list, "r") as f:
        line_counter = 0
        for line in f:
            tax_map[str(line_counter)] = line.strip()
            line_counter += 1
    print(tax_map)
    full_tree.rename_nodes(tax_map);
    if(hide_prefix):
        full_tree.write_tree_newick(output_file, hide_rooted_prefix=True)
    else:
        full_tree.write_tree_newick(output_file)

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--tax-list", required=True, type=click.Path(exists=True), help="Taxlist for the input tree")
@click.option("--output-file", required=True, type=str, help="Output file path for the relabeled tree")
@click.option("--hide-prefix", required=False, is_flag=True, help="Whether to include the rooted prefix in the tree file or not")
def relabel_tree(input_tree, tax_list, output_file, hide_prefix):
    """This script relabel the tree according to the input taxlist
    """
    relabel_tree_helper(input_tree, tax_list, output_file, hide_prefix)


if __name__ == "__main__":
    relabel_tree()
