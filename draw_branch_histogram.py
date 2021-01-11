import sys

import click
import dendropy
import matplotlib.pyplot as plt
import numpy as np

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="Input newick file containing the branch lengths.")
@click.option("--output-path", required=True, type=click.Path(exists=False), help="Path to output histogram plot of branch lengths.")
def draw_branch_histogram(input_tree, output_path):
    '''This program takes in a newick tree and writes a histogram of bootstrap support values for the nodes
    '''
    draw_branch_histograms_helper(input_tree, output_path)

def draw_branch_histogram_helper(input_tree, output_path, save_plot=True):
    current_tree = dendropy.Tree.get(path=input_tree, schema="newick")
    current_tree.encode_bipartitions()
    branch_lengths = []

    for edge in current_tree.postorder_edge_iter():
        if(edge.length is not None):
            if(float(edge.length) < 0.01):
                print(edge.length)
            branch_lengths.append(float(edge.length))
    bins = np.arange(0, 1.01, 0.01)
    if(save_plot):
        plt.clf()
        plt.hist(branch_lengths, bins=bins, density=False)
        plt.ylabel("Count")
        plt.xlabel("Branch Lengths");
        plt.savefig(output_path, bbox_inches='tight')

if __name__ == "__main__":
    draw_branch_histogram()
