import sys

import click
import dendropy
import matplotlib.pyplot as plt
import numpy as np

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="Input newick file containing the bootstrap support values.")
@click.option("--output-path", required=True, type=click.Path(exists=False), help="Path to output histogram plot of bootstrap support values.")
def draw_histograms(input_tree, output_path):
    '''This program takes in a newick tree and writes a histogram of bootstrap support values for the nodes
    '''
    draw_histograms_helper(input_tree, output_path)

def draw_histograms_helper(input_tree, output_path, save_plot=True):
    support_tree = dendropy.Tree.get(path=input_tree, schema="newick")
    support_tree.encode_bipartitions()
    support_values = []
    count_98 = 0
    count_99 = 0

    for node in support_tree.postorder_node_iter():
        if(node.label is not None):
            if(int(node.label) == 98):
                count_98 += 1
            elif(int(node.label) == 99):
                count_99 += 1
            support_values.append(float(node.label))
    bins = np.arange(0, 100, 1)
    if(save_plot):
        plt.clf()
        plt.hist(support_values, density=False, bins=bins, color="b")
        plt.ylabel("Count")
        plt.xlabel("Support Values");
        plt.savefig(output_path, bbox_inches='tight')
    return support_values


if __name__ == "__main__":
    draw_histograms()
