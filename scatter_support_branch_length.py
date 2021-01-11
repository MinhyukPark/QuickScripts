import sys

import click
import dendropy
import matplotlib.pyplot as plt
import numpy as np

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="Input newick file containing the bootstrap support values and branch lengths.")
@click.option("--output-path", required=True, type=click.Path(exists=False), help="Path to output scatter plot of bootstrap support values and branch lengths.")
def scatter_support_branch_length(input_tree, output_path):
    '''This program takes in a newick tree and writes a histogram of bootstrap support values for the nodes
    '''
    scatter_support_branch_length_helper(input_tree, output_path)

def scatter_support_branch_length_helper(input_tree, output_path, save_plot=True):
    current_tree = dendropy.Tree.get(path=input_tree, schema="newick")
    current_tree.encode_bipartitions()
    support_values = []
    branch_lengths = []

    for current_edge in current_tree.postorder_edge_iter():
        if(current_edge.head_node.taxon is None):
            if(current_edge.head_node.label is not None and current_edge.length is not None):
                support_values.append(float(current_edge.head_node.label))
                branch_lengths.append(float(current_edge.length))
    if(save_plot):
        plt.clf()
        plt.scatter(support_values, branch_lengths)
        plt.ylabel("Branch Lengths")
        plt.xlabel("Support Values");
        plt.ylim((0,1))
        plt.savefig(output_path, bbox_inches='tight')

if __name__ == "__main__":
    scatter_support_branch_length()
