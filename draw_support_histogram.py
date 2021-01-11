import sys

import click
import dendropy
import matplotlib.pyplot as plt
import numpy as np

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="Input newick file containing the bootstrap support values.")
@click.option("--output-path", required=True, type=click.Path(exists=False), help="Path to output histogram plot of bootstrap support values.")
def draw_support_histogram(input_tree, output_path):
    '''This program takes in a newick tree and writes a histogram of bootstrap support values for the nodes
    '''
    draw_histograms_helper(input_tree, output_path)

def draw_support_histogram_helper(input_tree, output_path, save_plot=True):
    support_tree = dendropy.Tree.get(path=input_tree, schema="newick")
    support_tree.encode_bipartitions()
    support_values = []

    for node in support_tree.postorder_node_iter():
        if(node.label is not None):
            support_values.append(float(node.label))
    # bins = np.arange(0, 101, 1)
    bins = np.arange(0, 1.01, 0.01)
    if(save_plot):
        plt.clf()
        plt.hist(support_values, density=False, bins=bins)
        plt.ylabel("Count")
        plt.xlabel("Support Values");
        plt.savefig(output_path, bbox_inches='tight')
    return support_values


if __name__ == "__main__":
    draw_support_histogram()
