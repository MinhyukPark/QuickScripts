import pathlib
import numpy as np

from pathlib import Path


from compare_trees import main as compare_trees


SCRATCH = "/u/sciteam/minhyuk2/scratch/"
SCRATCH_INPUT = SCRATCH + "input/"
CLUSTER_RUNS = "/u/sciteam/minhyuk2/cluster_runs/"

CLUSTER_RUN_MAP = {
    "nt78K": CLUSTER_RUNS + "GTM-experiment1/nt78K/INCOMPLETE/REPLICATE/",
    "100k_seeded": CLUSTER_RUNS + "GTM-experiment1/100k_seeded/INCOMPLETE/REPLICATE/",
    "1000M1-L": CLUSTER_RUNS + "GTM-experiment1/1000M1/Longest/INCOMPLETE/REPLICATE/",
    "1000M1-R": CLUSTER_RUNS + "GTM-experiment1/1000M1/random/INCOMPLETE/REPLICATE/",
}

MODEL_TREE_MAP = {
    "nt78K": SCRATCH_INPUT + "nt78K/huge.REPLICATE.sim.tree",
    "100k_seeded": SCRATCH_INPUT + "100k_seeded/REPLICATE/model/true.tt",
    "1000M1-L": SCRATCH_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
    "1000M1-R": SCRATCH_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
}

REPLICATE_MAP = {
    "nt78K": ["1", "2", "3", "4", "5"],
    "100k_seeded": ["0", "1", "2", "3", "4"],
    "1000M1-L": ["R0", "R1", "R2", "R3", "R4"],
    "1000M1-R": ["R0", "R1", "R2", "R3", "R4"],
}


def get_time_starting(cluster_run_path):
    starting_tree_file = cluster_run_path + "errors/fasttree.err"
    with open(starting_tree_file, "r") as f:
        for line in f:
            if "Total time:" in line:
                time_line = float(line.split(" ")[2])
                return time_line * (1/60) * (1/60) # hours


def get_time_guide(cluster_run_path):
    guide_tree_file = cluster_run_path + "output/iqtree2-incomplete.out"
    return get_time_iqtree(guide_tree_file)


def get_time_merge(cluster_run_path):
    gtm_file = cluster_run_path + "output/gtm.out"
    with open(gtm_file, "r") as f:
        for line in f:
            if "Finished GTM in" in line:
                time_line = float(line.split(" ")[3])
                print(line)
                print(time_line)
                return time_line * (1/60) * (1/60) # hours


def get_time_constraints(cluster_run_path):
    num_clusters = -42
    with open(cluster_run_path + "output/decompose.out", "r") as f:
        num_clusters = int(f.readlines()[0])
    current_max_time = 0
    for i in range(num_clusters):
        current_iqtree_path = cluster_run_path + "output/iqtree2-" + str(i) + ".out"
        current_max_time = max(current_max_time, get_time_iqtree(current_iqtree_path))
    return current_max_time


def get_time_iqtree(iqtree_path):
    with open(iqtree_path, "r") as f:
        for line in f:
            if "Total wall-clock time used:" in line:
                time_line = float(line.split(" ")[4])
                return time_line * (1/60) * (1/60) # hours



for dataset in CLUSTER_RUN_MAP:
    for incomplete in ["1", "5", "25"]:
        cumulative_fn_rate = 0.0
        time_starting_arr = []
        time_constraints_arr = []
        time_guide_arr = []
        time_merge_arr = []
        missing_count = 0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_cluster_run_path = CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("INCOMPLETE", incomplete)
            current_tree = current_cluster_run_path + "output/gtm.tree"

            if not Path(current_tree).is_file():
                missing_count += 1
                print(current_tree + " is missing")
                continue

            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
            if(fp != fn):
                print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                print("error: fp != fn")
            cumulative_fn_rate += (fn / ei1)
            time_starting_arr.append(get_time_starting(current_cluster_run_path))
            time_constraints_arr.append(get_time_constraints(current_cluster_run_path))
            time_guide_arr.append(get_time_guide(current_cluster_run_path))
            time_merge_arr.append(get_time_merge(current_cluster_run_path))
        if(missing_count == len(REPLICATE_MAP[dataset])):
            print(dataset + "-" + incomplete + " is not done yet.")
            print(str(missing_count) + " numbers of replicates are missing.")
        else:
            average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
            print(dataset + "-" + incomplete + ": " + str(average_fn_rate))
            print("starting tree time: " + str(np.median(time_strating_arr)))
            print("constraints tree time: " + str(np.median(time_constraints_arr)))
            print("guide tree time: " + str(np.median(time_guide_arr)))
            print("merge tree time: " + str(np.median(time_merge_arr)))

