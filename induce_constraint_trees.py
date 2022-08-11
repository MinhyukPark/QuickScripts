import random
import sys

import click
import treeswift
from Bio import SeqIO

from induce_tree import induce_tree_helper

PROJECTS = "/projects/tallis/minhyuk2/"
PROJECTS_INPUT = PROJECTS + "input/"
CLUSTER_RUNS = PROJECTS + "cluster_runs/"

EXPERIMENT_5_CLUSTER_RUN_MAP = {
    "nt10K": CLUSTER_RUNS + "GTM-experiment5/nt10K/Internal/REPLICATE/",
    "nt78K": CLUSTER_RUNS + "GTM-experiment5/nt78K/Internal/REPLICATE/",
    "100k_seeded": CLUSTER_RUNS + "GTM-experiment5/100k_seeded/Internal/REPLICATE/",
    "1000M1": CLUSTER_RUNS + "GTM-experiment5/1000M1/Internal/REPLICATE/",
}

EXPERIMENT_6_CLUSTER_RUN_MAP = {
    # "nt10K": CLUSTER_RUNS + "GTM-experiment6/nt10K/METHOD/REPLICATE/",
    # "nt78K": CLUSTER_RUNS + "GTM-experiment6/nt78K/METHOD/REPLICATE/",
    # "100k_seeded": CLUSTER_RUNS + "GTM-experiment6/100k_seeded/METHOD/REPLICATE/",
    "1000M1": CLUSTER_RUNS + "GTM-experiment6/1000M1/METHOD/REPLICATE/",
}

MODEL_TREE_MAP = {
    "nt10K": PROJECTS_INPUT + "nt10K/huge.REPLICATE.sim.tree",
    "nt78K": PROJECTS_INPUT + "nt78K/huge.REPLICATE.sim.single.tree",
    "100k_seeded": PROJECTS_INPUT + "100k_seeded/REPLICATE/model/true.tt",
    "1000M1-L": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
    "1000M1-R": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
    "1000M1": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
}

REPLICATE_MAP = {
    "nt10K": ["1", "2", "3", "4", "5"],
    "nt78K": ["1", "2", "3", "4", "5"],
    "100k_seeded": ["0", "1", "2", "3", "4"],
    "1000M1-L": ["R0", "R1", "R2", "R3", "R4"],
    "1000M1-R": ["R0", "R1", "R2", "R3", "R4"],
    "1000M1": ["R0", "R1", "R2", "R3", "R4"],
}

# @click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
# @click.option("--sequence-file", required=True, type=click.Path(exists=True), help="Sequence file with taxa labels to induce with")
# @click.option("--output-file", required=True, type=str, help="Output file path for the induced subtree")
# @click.option("--hide-prefix", required=False, is_flag=True, help="Whether to include the rooted prefix in the tree file or not")
@click.command()
def induce_constraint_trees():
    """This script induces the constraint trees
    """
    induce_1000M1_high_support()
    # for dataset in EXPERIMENT_5_CLUSTER_RUN_MAP:
    #     for replicate in REPLICATE_MAP[dataset]:
    #         cluster_run_path = EXPERIMENT_5_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate)
    #         current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
    #         num_clusters = -42
    #         with open(cluster_run_path + "output/decompose.out", "r") as f:
    #             num_clusters = int(f.readlines()[0])
    #         for i in range(num_clusters):
    #             current_sequence_file = cluster_run_path + "output/sequence_partition_" + str(i) + ".out"
    #             current_output_file = cluster_run_path + "output/sequence_partition_" + str(i) + ".model.tree"
    #             induce_tree_helper(current_model_tree, current_sequence_file, current_output_file, True)

    # for dataset in EXPERIMENT_6_CLUSTER_RUN_MAP:
    #     for replicate in REPLICATE_MAP[dataset]:
    #         for method in ["RAxML-ng", "FastTree"]:
    #             cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", method)
    #             current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
    #             num_clusters = -42
    #             with open(cluster_run_path + "output/decompose.out", "r") as f:
    #                 num_clusters = int(f.readlines()[0])
    #             for i in range(num_clusters):
    #                 current_sequence_file = cluster_run_path + "output/sequence_partition_" + str(i) + ".out"
    #                 current_output_file = cluster_run_path + "output/sequence_partition_" + str(i) + ".model.tree"
    #                 induce_tree_helper(current_model_tree, current_sequence_file, current_output_file, True)

def induce_1000M1_with_iqtree():
    for method in ["IQTree2"]:#["FastTree", "RAxML-ng", "IQTree2"]:
        for replicate in REPLICATE_MAP["1000M1"]:
            cluster_run_path = CLUSTER_RUNS + "GTM-experiment6/1000M1_with_IQTree_starting_tree/METHOD/REPLICATE/".replace("METHOD", method).replace("REPLICATE", replicate)
            current_model_tree = MODEL_TREE_MAP["1000M1"].replace("REPLICATE", replicate)
            num_clusters = 0
            print(cluster_run_path)
            with open(cluster_run_path + "output/decompose.out", "r") as f:
                num_clusters = int(f.readlines()[0])
            for i in range(num_clusters):
                current_sequence_file = cluster_run_path + "output/sequence_partition_" + str(i) + ".out"
                current_output_file = cluster_run_path + "output/sequence_partition_" + str(i) + ".model.tree"
                induce_tree_helper(current_model_tree, current_sequence_file, current_output_file, True)

def induce_1000M1_high_support():
    for method in ["FastTree", "RAxML-ng", "IQTree2"]:
        for replicate in REPLICATE_MAP["1000M1"]:
            cluster_run_path = CLUSTER_RUNS + "GTM-experiment6/1000M1_95_decompose/METHOD/REPLICATE/".replace("METHOD", method).replace("REPLICATE", replicate)
            current_model_tree = MODEL_TREE_MAP["1000M1"].replace("REPLICATE", replicate)
            num_clusters = 0
            print(cluster_run_path)
            with open(cluster_run_path + "output/decompose.out", "r") as f:
                num_clusters = int(f.readlines()[0])
            for i in range(num_clusters):
                current_sequence_file = cluster_run_path + "output/sequence_partition_" + str(i) + ".out"
                current_output_file = cluster_run_path + "output/sequence_partition_" + str(i) + ".model.tree"
                induce_tree_helper(current_model_tree, current_sequence_file, current_output_file, True)


if __name__ == "__main__":
    induce_constraint_trees()
