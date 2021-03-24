import sys

import click
import matplotlib.pyplot as plt
import numpy as np

@click.command()
@click.option("--input-dir", required=True, type=click.Path(exists=True), help="Input directory containing decomposed subsets.")
@click.option("--output-file", required=True, type=click.Path(exists=False), help="Path to output histogram plot.")
def get_cluster_size_histogram(input_dir, output_file):
    '''This program outputs a histogram of decomposed subset sizes
    '''
    subset_sizes = []

    num_subsets = int(open(input_dir + "decompose.out").readlines()[0])
    for i in range(num_subsets):
        current_subset_file = input_dir + "sequence_partition_" + str(i) + ".out"
        current_subset_size = 0
        with open(current_subset_file, "r") as f:
            for line in f:
                if line[0] == ">":
                    current_subset_size += 1
        subset_sizes.append(current_subset_size)


    bins = np.arange(0, 250, 1)
    plt.clf()
    plt.hist(subset_sizes, density=False, bins=bins, color="b")
    plt.ylabel("Count")
    plt.xlabel("Subset Sizes");
    plt.savefig(output_file, bbox_inches='tight')

if __name__ == "__main__":
    get_cluster_size_histogram()
