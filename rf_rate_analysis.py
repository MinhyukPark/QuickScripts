import sys
import pathlib
import numpy as np

from pathlib import Path


from compare_trees import main as compare_trees


PROJECTS = "/projects/tallis/minhyuk2/"
PROJECTS_INPUT = PROJECTS + "input/"
CLUSTER_RUNS = PROJECTS + "cluster_runs/"

EXPERIMENT_1_CLUSTER_RUN_MAP = {
    "nt10K": CLUSTER_RUNS + "GTM-experiment4/nt10K/INCOMPLETE/REPLICATE/",
}
    # "nt78K": CLUSTER_RUNS + "GTM-experiment1/nt78K/INCOMPLETE/REPLICATE/",
    # "100k_seeded": CLUSTER_RUNS + "GTM-experiment1/100k_seeded/INCOMPLETE/REPLICATE/",
    # "1000M1-L": CLUSTER_RUNS + "GTM-experiment1/1000M1/Longest/INCOMPLETE/REPLICATE/",
    # "1000M1-R": CLUSTER_RUNS + "GTM-experiment1/1000M1/random/INCOMPLETE/REPLICATE/",
# }

EXPERIMENT_2_CLUSTER_RUN_MAP = {
    "nt10K": CLUSTER_RUNS + "GTM-experiment4/nt10K/GUIDETREE/REPLICATE/",
}
    # "nt78K": CLUSTER_RUNS + "GTM-experiment2/nt78K/GUIDETREE/REPLICATE/",
    # "100k_seeded": CLUSTER_RUNS + "GTM-experiment2/100k_seeded/GUIDETREE/REPLICATE/",
    # "1000M1": CLUSTER_RUNS + "GTM-experiment2/1000M1/GUIDETREE/REPLICATE/",
# }

EXPERIMENT_4_CLUSTER_RUN_MAP = {
    "nt10K": CLUSTER_RUNS + "GTM-experiment4/nt10K/SUBSET/REPLICATE/",
    "nt78K": CLUSTER_RUNS + "GTM-experiment4/nt78K/SUBSET/REPLICATE/",
    "100k_seeded": CLUSTER_RUNS + "GTM-experiment4/100k_seeded/SUBSET/REPLICATE/",
}

MODEL_TREE_MAP = {
    "nt78K": PROJECTS_INPUT + "nt78K/huge.REPLICATE.sim.tree",
    "100k_seeded": PROJECTS_INPUT + "100k_seeded/REPLICATE/model/true.tt",
    "1000M1-L": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
    "1000M1-R": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
    "1000M1": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
}

REPLICATE_MAP = {
    "nt78K": ["1", "2", "3", "4", "5"],
    "100k_seeded": ["0", "1", "2", "3", "4"],
    "1000M1-L": ["R0", "R1", "R2", "R3", "R4"],
    "1000M1-R": ["R0", "R1", "R2", "R3", "R4"],
    "1000M1": ["R0", "R1", "R2", "R3", "R4"],
}


def get_time_starting(cluster_run_path):
    starting_tree_file = cluster_run_path + "errors/fasttree.err"
    return get_time_fasttree(starting_tree_file)

def get_time_fasttree(starting_tree_file):
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

def get_gtm_experiment_1():
    for dataset in EXPERIMENT_1_CLUSTER_RUN_MAP:
        for incomplete in ["1", "5", "25"]:
            sys.stdout.flush()
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0
            cumulative_starting_fn_rate = 0.0
            cumulative_starting_fp_rate = 0.0
            time_starting_arr = []
            time_constraints_arr = []
            time_guide_arr = []
            time_merge_arr = []
            missing_count = 0
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_cluster_run_path = EXPERIMENT_1_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("INCOMPLETE", incomplete)
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
                cumulative_fp_rate += (fp / ei1)
                time_starting_arr.append(get_time_starting(current_cluster_run_path))
                time_constraints_arr.append(get_time_constraints(current_cluster_run_path))
                time_guide_arr.append(get_time_guide(current_cluster_run_path))
                time_merge_arr.append(get_time_merge(current_cluster_run_path))
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + incomplete + " is not done yet.")
                print(str(missing_count) + " numbers of replicates are missing.")
            else:
                average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                print(dataset + "-" + incomplete + " gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + incomplete + " gtm fp rate: " + str(average_fp_rate))
                print("starting tree time: " + str(np.median(time_starting_arr)))
                print("constraints tree time: " + str(np.median(time_constraints_arr)))
                print("guide tree time: " + str(np.median(time_guide_arr)))
                print("merge tree time: " + str(np.median(time_merge_arr)))

def get_gtm_experiment_2():
    for dataset in EXPERIMENT_2_CLUSTER_RUN_MAP:
        sys.stdout.flush()
        cumulative_fn_rate = 0
        cumulative_fp_rate = 0
        time_starting_arr = []
        time_constraints_arr = []
        time_merge_arr = []
        missing_count = 0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_cluster_run_path = EXPERIMENT_2_CLUSTER_RUN_MAP[dataset].replace("GUIDETREE", "FastTree").replace("REPLICATE", replicate)
            current_output_tree = current_cluster_run_path + "output/gtm.tree"
            if not Path(current_output_tree).is_file():
                missing_count += 1
                print(current_output_tree + " is missing")
                continue
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_output_tree)
            if(fp != fn):
                print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                print("error: fp != fn")
            cumulative_fn_rate += (fn / ei1)
            cumulative_fp_rate += (fp / ei1)
            time_starting_arr.append(get_time_starting(current_cluster_run_path))
            time_constraints_arr.append(get_time_constraints(current_cluster_run_path))
            time_merge_arr.append(get_time_merge(current_cluster_run_path))
        if(missing_count == len(REPLICATE_MAP[dataset])):
            print(dataset + "-" + incomplete + " is not done yet.")
            print(str(missing_count) + " numbers of replicates are missing.")
        else:
            average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
            average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
            print(dataset + " gtm fn rate: " + str(average_fn_rate))
            print(dataset + " gtm fp rate: " + str(average_fp_rate))
            print("starting tree time: " + str(np.median(time_starting_arr)))
            print("constraints tree time: " + str(np.median(time_constraints_arr)))
            print("guide tree time is the same as the starting tree time")
            print("merge tree time: " + str(np.median(time_merge_arr)))

        time_starting_arr = []
        time_constraints_arr = []
        time_guide_arr = []
        time_merge_arr = []
        missing_count = 0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_cluster_run_path = EXPERIMENT_2_CLUSTER_RUN_MAP[dataset].replace("GUIDETREE", "Internal").replace("REPLICATE", replicate)
            current_output_tree = current_cluster_run_path + "output/gtm.tree"
            if not Path(current_output_tree).is_file():
                missing_count += 1
                print(current_output_tree + " is missing")
                continue
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_output_tree)
            if(fp != fn):
                print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                print("error: fp != fn")
            cumulative_fn_rate += (fn / ei1)
            cumulative_fp_rate += (fp / ei1)
            time_starting_arr.append(get_time_starting(current_cluster_run_path))
            time_constraints_arr.append(get_time_constraints(current_cluster_run_path))
            time_merge_arr.append(get_time_merge(current_cluster_run_path))
            time_guide_arr.append(get_time_guide(current_cluster_run_path))
        if(missing_count == len(REPLICATE_MAP[dataset])):
            print(dataset + "-" + incomplete + " is not done yet.")
            print(str(missing_count) + " numbers of replicates are missing.")
        else:
            average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
            average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
            print(dataset + "-" + incomplete + " gtm fn rate: " + str(average_fn_rate))
            print(dataset + "-" + incomplete + " gtm fp rate: " + str(average_fp_rate))
            print("starting tree time: " + str(np.median(time_starting_arr)))
            print("constraints tree time: " + str(np.median(time_constraints_arr)))
            print("guide tree time: " + str(np.median(time_guide_arr)))
            print("merge tree time: " + str(np.median(time_merge_arr)))

def get_fasttree_and_internal():
    for dataset in EXPERIMENT_2_CLUSTER_RUN_MAP:
        cumulative_fasttree_fn_rate = 0.0
        cumulative_fasttree_fp_rate = 0.0
        time_fasttree_arr = []
        fasttree_missing_count = 0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_cluster_run_path = EXPERIMENT_2_CLUSTER_RUN_MAP[dataset].replace("GUIDETREE", "FastTree").replace("REPLICATE", replicate)
            current_fasttree_tree = current_cluster_run_path + "output/fasttree.out"
            if not Path(current_fasttree_tree).is_file():
                fasttree_missing_count += 1
                print(current_fasttree_tree + " is missing")
                continue
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_fasttree_tree)
            if(fp != fn):
                print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                print("error: fp != fn")
            cumulative_fasttree_fn_rate += (fn / ei1)
            cumulative_fasttree_fp_rate += (fp / ei1)
            time_fasttree_arr.append(get_time_fasttree(current_cluster_run_path + "errors/fasttree.err"))

        average_fasttree_fn_rate = cumulative_fasttree_fn_rate / (len(REPLICATE_MAP[dataset]) - fasttree_missing_count)
        average_fasttree_fp_rate = cumulative_fasttree_fp_rate / (len(REPLICATE_MAP[dataset]) - fasttree_missing_count)
        print(dataset + "-fasttree fn rate: " + str(average_fasttree_fn_rate))
        print(dataset + "-fasttree fp rate: " + str(average_fasttree_fp_rate))
        print("fasttree tree time: " + str(np.median(time_fasttree_arr)))

        cumulative_internal_fn_rate = 0.0
        cumulative_internal_fp_rate = 0.0
        time_internal_arr = []
        internal_missing_count = 0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_cluster_run_path = EXPERIMENT_2_CLUSTER_RUN_MAP[dataset].replace("GUIDETREE", "Internal").replace("REPLICATE", replicate)
            current_internal_tree = current_cluster_run_path + "output/fasttree-int.out"
            if not Path(current_internal_tree).is_file():
                internal_missing_count += 1
                print(current_internal_tree + " is missing")
                continue
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_internal_tree)
            if(fp != fn):
                print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                print("error: fp != fn")
            cumulative_internal_fn_rate += (fn / ei1)
            cumulative_internal_fp_rate += (fp / ei1)
            time_internal_arr.append(get_time_fasttree(current_cluster_run_path + "errors/fasttree-int.err"))

        average_internal_fn_rate = cumulative_internal_fn_rate / (len(REPLICATE_MAP[dataset]) - internal_missing_count)
        average_internal_fp_rate = cumulative_internal_fp_rate / (len(REPLICATE_MAP[dataset]) - internal_missing_count)
        print(dataset + "-internal fn rate: " + str(average_internal_fn_rate))
        print(dataset + "-internal fp rate: " + str(average_internal_fp_rate))
        print("internal tree time: " + str(np.median(time_internal_arr)))

def get_incomplete_tree_error():
    for dataset in EXPERIMENT_1_CLUSTER_RUN_MAP:
        for incomplete in ["1", "5", "25"]:
            sys.stdout.flush()
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0
            missing_count = 0
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_cluster_run_path = EXPERIMENT_1_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("INCOMPLETE", incomplete)
                current_induced_model_tree = current_cluster_run_path + "output/incomplete-induced-model.tree"
                current_incomplete_tree = current_cluster_run_path + "output/guide.treefile"

                if not Path(current_incomplete_tree).is_file():
                    missing_count += 1
                    print(current_incomplete_tree + " is missing")
                    continue

                nl, ei1, ei2, fp, fn, rf = compare_trees(current_induced_model_tree, current_incomplete_tree)
                if(fp != fn):
                    print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    print("error: fp != fn")
                cumulative_fn_rate += (fn / ei1)
                cumulative_fp_rate += (fp / ei1)
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + incomplete + " is not done yet.")
                print(str(missing_count) + " numbers of replicates are missing.")
            else:
                average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                print(dataset + "-" + incomplete + " gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + incomplete + " gtm fp rate: " + str(average_fp_rate))


def get_gtm_experiment_4():
    for dataset in EXPERIMENT_4_CLUSTER_RUN_MAP:
        for subset in ["250", "1000"]:
            sys.stdout.flush()
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0
            cumulative_starting_fn_rate = 0.0
            cumulative_starting_fp_rate = 0.0
            time_starting_arr = []
            time_constraints_arr = []
            time_guide_arr = []
            time_merge_arr = []
            missing_count = 0
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_cluster_run_path = EXPERIMENT_4_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("SUBSET", subset)
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
                cumulative_fp_rate += (fp / ei1)
                time_starting_arr.append(get_time_starting(current_cluster_run_path))
                time_constraints_arr.append(get_time_constraints(current_cluster_run_path))
                time_guide_arr.append(get_time_guide(current_cluster_run_path))
                time_merge_arr.append(get_time_merge(current_cluster_run_path))
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + incomplete + " is not done yet.")
                print(str(missing_count) + " numbers of replicates are missing.")
            else:
                average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                print(dataset + "-" + incomplete + " gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + incomplete + " gtm fp rate: " + str(average_fp_rate))
                print("starting tree time: " + str(np.median(time_starting_arr)))
                print("constraints tree time: " + str(np.median(time_constraints_arr)))
                print("guide tree time: " + str(np.median(time_guide_arr)))
                print("merge tree time: " + str(np.median(time_merge_arr)))


# get_incomplete_tree_error()
print("EXPERIMENT 1")
# get_fasttree_and_internal()
get_gtm_experiment_1()
print("EXPERIMENT 2")
get_gtm_experiment_2()
print("EXPERIMENT 4")
get_gtm_experiment_4()
