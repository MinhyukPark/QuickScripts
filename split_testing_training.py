import random
import sys

import click
import treeswift
from Bio import SeqIO

@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--sequence-file", required=True, type=click.Path(exists=True), help="Aligned sequence file on the full taxa")
@click.option("--output-prefix", required=True, type=str, help="Output file prefix for each subset")
@click.option("--minimum-size", required=False, type=int, help="Minimum size of each set")
def split_testing_training(input_tree, sequence_file, output_prefix, minimum_size):
    """This script decomposes the input tree and returns two clades by splitting at two random edges
    """
    f = open("/projects/tallis/minhyuk2/input/PASTA/RNASim/1000000/p.out", "a+")
    f.write("starting split_testing_training\n")
    f.close()
    guide_tree = treeswift.read_tree_newick(input_tree)
    num_nodes = guide_tree.num_nodes()
    testing_obtained = False
    training_obtained = False
    f = open("/projects/tallis/minhyuk2/input/PASTA/RNASim/1000000/p.out", "a+")
    f.write("starting split_testing_training\n")
    f.close()
    while(not testing_obtained or not training_obtained):
        testing_obtained = False
        training_obtained = False
        training_subtree_index = random.randrange(num_nodes / 2)
        testing_subtree_index = random.randrange(num_nodes / 2, num_nodes)
        for node_index,node in enumerate(tree.traverse_postorder(leaves=False, internal=True)):
            if(node_index == training_subtree_index):
                training_tree = extract_subtree(node)
                if(training_tree.num_num_nodes() < minimum_size):
                    break
                else:
                    training_obtained = True
            elif(node_index == testing_subtree_index):
                testing_tree = extract_subtree(node)
                if(testing_tree.num_num_nodes() < minimum_size):
                    break
                else:
                    testing_obtained = True

    f = open("/projects/tallis/minhyuk2/input/PASTA/RNASim/1000000/p.out", "a+")
    f.write("obtained both trees. partitioning.\n")
    f.close()

    current_sequence_partition = list(testing_tree.traverse_postorder(leaves=True, internal=False))
    for sequence in SeqIO.parse(open(sequence_file), "fasta"):
        if(sequence.id in current_cluster):
            current_sequence_partition.append(sequence)
    SeqIO.write(current_sequence_partition, output_prefix + "testing.out", "fasta")

    f = open("/projects/tallis/minhyuk2/input/PASTA/RNASim/1000000/p.out", "a+")
    f.write("partitioned testing\n")
    f.close()

    current_sequence_partition = list(training_tree.traverse_postorder(leaves=True, internal=False))
    for sequence in SeqIO.parse(open(sequence_file), "fasta"):
        if(sequence.id in current_cluster):
            current_sequence_partition.append(sequence)
    SeqIO.write(current_sequence_partition, output_prefix + "training.out", "fasta")

    f = open("/projects/tallis/minhyuk2/input/PASTA/RNASim/1000000/p.out", "a+")
    f.write("partitioned training\n")
    f.close()

if __name__ == "__main__":
    split_testing_training()
