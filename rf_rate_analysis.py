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
    # "nt10K": CLUSTER_RUNS + "GTM-experiment4/nt10K/SUBSET/REPLICATE/",
    # "nt78K": CLUSTER_RUNS + "GTM-experiment4/nt78K/SUBSET/REPLICATE/",
    # "100k_seeded": CLUSTER_RUNS + "GTM-experiment4/100k_seeded/SUBSET/REPLICATE/",
    "1000M1": CLUSTER_RUNS + "GTM-experiment4/1000M1/SUBSET/REPLICATE/",
}

EXPERIMENT_5_CLUSTER_RUN_MAP = {
    # "nt10K": CLUSTER_RUNS + "GTM-experiment5/nt10K/Internal/REPLICATE/",
    # "nt78K": CLUSTER_RUNS + "GTM-experiment5/nt78K/Internal/REPLICATE/",
    # "100k_seeded": CLUSTER_RUNS + "GTM-experiment5/100k_seeded/Internal/REPLICATE/",
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
    "nt78K": PROJECTS_INPUT + "nt78K/huge.REPLICATE.sim.tree",
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

TREEMERGE_1000M1_MAP = {
    "raxmlng": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment6/1000M1_with_TreeMerge/RAxML-ng/REPLICATE/output/",
    "iqtree2": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment6/1000M1_with_TreeMerge/IQTree2/REPLICATE/output/",
    "fasttree": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment6/1000M1_with_TreeMerge/FastTree/REPLICATE/output/",
}

BASE_1000M1_MAP = {
    # "raxmlng": "/projects/tallis/minhyuk2/cluster_runs/individual_test/1000M1RAxML/REPLICATE/raxmlng.raxml.treefileTODO",
    "raxmlng": "/projects/tallis/minhyuk2/cluster_runs/individual_test/1000M1RAxML/REPLICATE/output.raxml.lastTree.TMP",
    # "fasttree": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment4/1000M1/250/REPLICATE/output/fasttree.out",
    # "iqtree2": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment4/1000M1/250/REPLICATE/output/guide.treefile",
}

BASE_1000M1_ERR_MAP = {
    # "raxmlng": "/projects/tallis/minhyuk2/cluster_runs/individual_test/1000M1RAxML/REPLICATE/raxmlng.raxml.log",
    # "fasttree": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment4/1000M1/250/REPLICATE/errors/fasttree.err",
    # "iqtree2": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment4/1000M1/250/REPLICATE/output/guide.log",
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


def get_time_constraints(cluster_run_path, method="IQTree2"):
    num_clusters = -42
    with open(cluster_run_path + "output/decompose.out", "r") as f:
        num_clusters = int(f.readlines()[0])
    current_max_time = 0
    if(method == "IQTree2"):
        for i in range(num_clusters):
            current_iqtree_path = cluster_run_path + "output/iqtree2-" + str(i) + ".out"
            current_max_time = max(current_max_time, get_time_iqtree(current_iqtree_path))
    elif(method == "FastTree"):
        for i in range(num_clusters):
            current_fasttree_path = cluster_run_path + "errors/fasttree-" + str(i) + ".err"
            current_max_time = max(current_max_time, get_time_fasttree(current_fasttree_path))
    elif(method == "RAxML-ng"):
        for i in range(num_clusters):
            current_raxmlng_path = cluster_run_path + "output/sequence_partition_" + str(i) + ".raxml.log"
            current_max_time = max(current_max_time, get_time_raxmlng(current_raxmlng_path))
    return current_max_time


def get_constraint_accuracies(cluster_run_path, constraint_method):
    num_clusters = -42
    with open(cluster_run_path + "output/decompose.out", "r") as f:
        num_clusters = int(f.readlines()[0])
    fn_rate = 0.0
    fp_rate = 0.0
    for i in range(num_clusters):
        current_model_tree = cluster_run_path + "output/sequence_partition_" + str(i) + ".model.tree"
        if(constraint_method == "FastTree"):
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, cluster_run_path + "output/fasttree-" + str(i) + ".out")
            fn_rate += (fn / ei1)
            fp_rate += (fp / ei1)
            print(f'fn_rate: {fn / ei1} and fp_rate: {fp / ei1}')
        elif(constraint_method == "RAxML-ng"):
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, cluster_run_path + "output/sequence_partition_" + str(i) + ".raxml.bestTree")
            fn_rate += (fn / ei1)
            fp_rate += (fp / ei1)
            print(f'fn_rate: {fn / ei1} and fp_rate: {fp / ei1}')
        elif(constraint_method == "IQTree2"):
            print(current_model_tree)
            print(cluster_run_path + "output/sequence_partition_" + str(i) + ".treefile")
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, cluster_run_path + "output/sequence_partition_" + str(i) + ".treefile")
            fn_rate += (fn / ei1)
            fp_rate += (fp / ei1)
            print(f'fn_rate: {fn / ei1} and fp_rate: {fp / ei1}')
        # if(fp != fn):
        #     print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
        #     print("error: fp != fn")

    fn_rate /= num_clusters
    fp_rate /= num_clusters
    return {
        "fn": fn_rate,
        "fp": fp_rate,
    }


def get_time_iqtree(iqtree_path):
    with open(iqtree_path, "r") as f:
        for line in f:
            if "Total wall-clock time used:" in line:
                time_line = float(line.split(" ")[4])
                return time_line * (1/60) * (1/60) # hours

def get_time_raxmlng(raxmlng_path):
    with open(raxmlng_path, "r") as f:
        for line in f:
            if "Elapsed time" in line:
                time_line = float(line.split(" ")[2])
                return time_line * (1/60) * (1/60)

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
        for subset in ["250", "500"]:
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
                if(dataset == "1000M1"):
                    time_guide_arr.append(get_time_iqtree(current_cluster_run_path + "output/iqtree2-complete.out"))
                else:
                    time_guide_arr.append(get_time_guide(current_cluster_run_path))
                time_merge_arr.append(get_time_merge(current_cluster_run_path))
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + subset + " is not done yet.")
                print(str(missing_count) + " numbers of replicates are missing.")
            else:
                average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                print(dataset + "-" + subset + "-gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + subset + "-gtm fp rate: " + str(average_fp_rate))
                print("starting tree time: " + str(np.median(time_starting_arr)))
                print("constraints tree time: " + str(np.median(time_constraints_arr)))
                print("guide tree time: " + str(np.median(time_guide_arr)))
                print("merge tree time: " + str(np.median(time_merge_arr)))

def get_gtm_experiment_5():
    for dataset in EXPERIMENT_5_CLUSTER_RUN_MAP:
        sys.stdout.flush()
        cumulative_fn_rate = 0.0
        cumulative_fp_rate = 0.0
        cumulative_starting_fn_rate = 0.0
        cumulative_starting_fp_rate = 0.0
        cumulative_constraint_fn_rate = 0.0
        cumulative_constraint_fp_rate = 0.0
        time_starting_arr = []
        time_constraints_arr = []
        time_guide_arr = []
        time_merge_arr = []
        missing_count = 0
        fn = 0.0
        fp = 0.0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_cluster_run_path = EXPERIMENT_5_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate)
            current_tree = current_cluster_run_path + "output/gtm.tree"

            if not Path(current_tree).is_file():
                missing_count += 1
                print(current_tree + " is missing")
                continue

            # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
            # if(fp != fn):
                # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                # print("error: fp != fn")
            # cumulative_fn_rate += (fn / ei1)
            # cumulative_fp_rate += (fp / ei1)

            print("dataset: " + dataset + " method: IQTree2 replicate: " + replicate)
            constraint_accuracies = get_constraint_accuracies(current_cluster_run_path, "IQTree2")
            # cumulative_constraint_fn_rate += constraint_accuracies["fn"]
            # cumulative_constraint_fp_rate += constraint_accuracies["fp"]

            # time_starting_arr.append(get_time_fasttree(current_cluster_run_path + "errors/fasttree-int.err"))
            # time_constraints_arr.append(get_time_constraints(current_cluster_run_path))
            # if(dataset == "1000M1"):
                # time_guide_arr.append(get_time_iqtree(current_cluster_run_path + "output/iqtree2-complete.out"))
            # else:
                # time_guide_arr.append(get_time_guide(current_cluster_run_path))
            # time_merge_arr.append(get_time_merge(current_cluster_run_path))
        if(missing_count == len(REPLICATE_MAP[dataset])):
            print(dataset + " is not done yet.")
            print(str(missing_count) + " numbers of replicates are missing.")
        else:
            print()
            # average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
            # average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
            # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
            # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
            # print(dataset + "-internal gtm fn rate: " + str(average_fn_rate))
            # print(dataset + "-internal gtm fp rate: " + str(average_fp_rate))
            # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
            # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
            # print("starting tree time: " + str(np.median(time_starting_arr)))
            # print("constraints tree time: " + str(np.median(time_constraints_arr)))
            # print("guide tree time: " + str(np.median(time_guide_arr)))
            # print("merge tree time: " + str(np.median(time_merge_arr)))

def get_1000M1_base():
    dataset = "1000M1"
    # for method in ["fasttree", "iqtree2"]:
    for method in ["raxmlng"]:
        cumulative_fn_rate = 0.0
        cumulative_fp_rate = 0.0
        cumulative_time = 0.0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_tree = BASE_1000M1_MAP[method].replace("REPLICATE", replicate)
            # error_path = BASE_1000M1_ERR_MAP[method].replace("REPLICATE", replicate)
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
            if(fp != fn):
                print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                print("error: fp != fn")
            cumulative_fn_rate += (fn / ei1)
            cumulative_fp_rate += (fp / ei1)
            if(method == "raxmlng"):
                # cumulative_time += get_time_raxmlng(error_path)
                pass
            elif(method == "iqtree2"):
                cumulative_time += get_time_iqtree(error_path)
            elif(method == "fasttree"):
                cumulative_time += get_time_fasttree(error_path)

        print(dataset + "-" + method + " fn rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
        print(dataset + "-" + method + " fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
        print(dataset + "-" + method + " average time: ", str(cumulative_time / len(REPLICATE_MAP[dataset])))


def get_gtm_experiment_6():
    for dataset in EXPERIMENT_6_CLUSTER_RUN_MAP:
        for constraint_method in ["FastTree", "RAxML-ng"]:
            sys.stdout.flush()
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0

            cumulative_constraint_fn_rate = 0.0
            cumulative_constraint_fp_rate = 0.0
            cumulative_guide_fn_rate = 0.0
            cumulative_guide_fp_rate = 0.0

            time_starting_arr = []
            time_constraints_arr = []
            time_guide_arr = []
            time_merge_arr = []
            missing_count = 0
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method)
                current_tree = current_cluster_run_path + "output/gtm.tree"

                if not Path(current_tree).is_file():
                    missing_count += 1
                    print(current_tree + " is missing")
                    continue

                # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                # if(fp != fn):
                #     print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                #     print("error: fp != fn")
                # cumulative_fn_rate += (fn / ei1)
                # cumulative_fp_rate += (fp / ei1)

                # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                # if(fp != fn):
                #     print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                #     print("error: fp != fn")
                # cumulative_guide_fn_rate += (fn / ei1)
                # cumulative_guide_fp_rate += (fp / ei1)

                print("dataset: " + dataset + " method: " + constraint_method  + " replicate: " + replicate)
                constraint_accuracies = get_constraint_accuracies(current_cluster_run_path, constraint_method)
                # cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                # cumulative_constraint_fp_rate += constraint_accuracies["fp"]

                # time_starting_arr.append(get_time_starting(current_cluster_run_path))
                # time_constraints_arr.append(get_time_constraints(current_cluster_run_path, method=constraint_method))
                # if(dataset == "1000M1"):
                    # time_guide_arr.append(get_time_iqtree(current_cluster_run_path + "output/iqtree2-complete.out"))
                # else:
                    # time_guide_arr.append(get_time_guide(current_cluster_run_path))
                # time_merge_arr.append(get_time_merge(current_cluster_run_path))
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + constraint_method + " is not done yet.")
                print(str(missing_count) + " numbers of replicates are missing.")
            else:
                print()
                # average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)

                # print(dataset + "-" + constraint_method + "-gtm fn rate: " + str(average_fn_rate))
                # print(dataset + "-" + constraint_method + "-gtm fp rate: " + str(average_fp_rate))
                # print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                # print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                # print("starting tree time: " + str(np.median(time_starting_arr)))
                # print("constraints tree time: " + str(np.median(time_constraints_arr)))
                # # print("guide tree time: " + str(np.median(time_guide_arr)))
                # print("merge tree time: " + str(np.median(time_merge_arr)))


def get_200M1_results():
    for method in ["IQTree2", "FastTree", "RAxML-ng"]:
        cumulative_fn_rate = 0.0
        cumulative_fp_rate = 0.0
        for replicate in ["R0", "R1", "R2", "R3", "R4"]:
            current_path = "/projects/tallis/minhyuk2/cluster_runs/individual_test/200M1Tests/" + replicate + "/"
            current_model_tree = current_path + "true_tree.tre"
            current_tree = None
            if(method == "IQTree2"):
                current_tree = current_path + "iqtree2output.treefile"
            elif(method == "FastTree"):
                current_tree = current_path + "fasttreeoutput.out"
            elif(method == "RAxML-ng"):
                current_tree = current_path + "raxmloutput.raxml.bestTree"

            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
            # if(fp != fn):
            #     print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
            #     print("error: fp != fn")
            cumulative_fn_rate += (fn / ei1)
            cumulative_fp_rate += (fp / ei1)

        average_fn_rate = cumulative_fn_rate / 5
        average_fp_rate = cumulative_fp_rate / 5
        print(method + " fn rate: " + str(average_fn_rate))
        print(method + " fp rate: " + str(average_fp_rate))


def get_1000M1_treemerge_results():
    dataset = "1000M1"
    for method in ["fasttree", "iqtree2", "raxmlng"]:
        for distance_type in ["mldist", "nodedist"]:
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0
            cumulative_time = 0.0
            missing_count = 0
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_tree = TREEMERGE_1000M1_MAP[method].replace("REPLICATE", replicate) + "treemerge_" + distance_type + ".tree"

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
                # if(method == "raxmlng"):
                #     cumulative_time += get_time_raxmlng(error_path)
                #     pass
                # elif(method == "iqtree2"):
                #     cumulative_time += get_time_iqtree(error_path)
                # elif(method == "fasttree"):
                #     cumulative_time += get_time_fasttree(error_path)
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + method + "-" + distance_type + " not done yet")
            else:
                print(dataset + "-" + method + "-" + distance_type + " fn rate: ", str(cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)))
                print(dataset + "-" + method + "-" + distance_type + " fp rate: ", str(cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)))
                print(dataset + "-" + method + "-" + distance_type + " average time: ", str(cumulative_time / (len(REPLICATE_MAP[dataset]) - missing_count)))

def get_gtm_experiment_6_iqtree_starting():
    for dataset in EXPERIMENT_6_CLUSTER_RUN_MAP:
        for constraint_method in ["IQTree2"]:#["FastTree", "RAxML-ng"]:
            sys.stdout.flush()
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0

            cumulative_constraint_fn_rate = 0.0
            cumulative_constraint_fp_rate = 0.0
            cumulative_guide_fn_rate = 0.0
            cumulative_guide_fp_rate = 0.0

            time_starting_arr = []
            time_constraints_arr = []
            time_guide_arr = []
            time_merge_arr = []
            missing_count = 0
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method).replace("1000M1", "1000M1_with_IQTree_starting_tree")
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

                nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                if(fp != fn):
                    print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    print("error: fp != fn")
                cumulative_guide_fn_rate += (fn / ei1)
                cumulative_guide_fp_rate += (fp / ei1)

                print("dataset: " + dataset + " method: " + constraint_method  + " replicate: " + replicate)
                sys.stdout.flush()
                constraint_accuracies = get_constraint_accuracies(current_cluster_run_path, constraint_method)
                cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                cumulative_constraint_fp_rate += constraint_accuracies["fp"]

                # time_starting_arr.append(get_time_iqtree(current_cluster_run_path))
                time_constraints_arr.append(get_time_constraints(current_cluster_run_path, method=constraint_method))
                if(dataset == "1000M1"):
                    time_guide_arr.append(get_time_iqtree(current_cluster_run_path + "output/iqtree2-complete.out"))
                else:
                    time_guide_arr.append(get_time_guide(current_cluster_run_path))
                time_merge_arr.append(get_time_merge(current_cluster_run_path))
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + constraint_method + " is not done yet.")
                print(str(missing_count) + " numbers of replicates are missing.")
            else:
                average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)

                print(dataset + "-" + constraint_method + "-gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + constraint_method + "-gtm fp rate: " + str(average_fp_rate))
                print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                # print("starting tree time: " + str(np.median(time_starting_arr)))
                print("constraints tree time: " + str(np.median(time_constraints_arr)))
                print("guide tree time: " + str(np.median(time_guide_arr)))
                print("merge tree time: " + str(np.median(time_merge_arr)))


def get_gtm_experiment_6_high_support():
    for dataset in EXPERIMENT_6_CLUSTER_RUN_MAP:
        for constraint_method in ["IQTree2", "FastTree", "RAxML-ng"]:
            sys.stdout.flush()
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0

            cumulative_constraint_fn_rate = 0.0
            cumulative_constraint_fp_rate = 0.0
            cumulative_guide_fn_rate = 0.0
            cumulative_guide_fp_rate = 0.0

            time_starting_arr = []
            time_constraints_arr = []
            time_guide_arr = []
            time_merge_arr = []
            missing_count = 0
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method).replace("1000M1", "1000M1_95_decompose")
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

                nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                if(fp != fn):
                    print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    print("error: fp != fn")
                cumulative_guide_fn_rate += (fn / ei1)
                cumulative_guide_fp_rate += (fp / ei1)

                print("dataset: " + dataset + " method: " + constraint_method  + " replicate: " + replicate)
                sys.stdout.flush()
                constraint_accuracies = get_constraint_accuracies(current_cluster_run_path, constraint_method)
                cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                cumulative_constraint_fp_rate += constraint_accuracies["fp"]

                time_starting_arr.append(get_time_starting(current_cluster_run_path))
                time_constraints_arr.append(get_time_constraints(current_cluster_run_path, method=constraint_method))
                if(dataset == "1000M1"):
                    time_guide_arr.append(get_time_iqtree(current_cluster_run_path + "output/iqtree2-complete.out"))
                else:
                    time_guide_arr.append(get_time_guide(current_cluster_run_path))
                time_merge_arr.append(get_time_merge(current_cluster_run_path))
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + constraint_method + " is not done yet.")
                print(str(missing_count) + " numbers of replicates are missing.")
            else:
                average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)

                print(dataset + "-" + constraint_method + "-gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + constraint_method + "-gtm fp rate: " + str(average_fp_rate))
                print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                print("starting tree time: " + str(np.median(time_starting_arr)))
                print("constraints tree time: " + str(np.median(time_constraints_arr)))
                print("guide tree time: " + str(np.median(time_guide_arr)))
                print("merge tree time: " + str(np.median(time_merge_arr)))



# get_incomplete_tree_error()
#print("EXPERIMENT 1")
# get_fasttree_and_internal()
# get_gtm_experiment_1()
#print("EXPERIMENT 2")
# get_gtm_experiment_2()
# print("EXPERIMENT 4")
# get_gtm_experiment_4()
# print("EXPERIMENT 5")
# get_gtm_experiment_5()
# print("GET 1000M1 Base")
# get_1000M1_base()
#get_200M1_results()
# print("1000M1 PARTIAL TREEMERGE RESULTS")
# get_1000M1_treemerge_results()

# print("EXPERIMENT 6 IQTREE STARTING")
# get_gtm_experiment_6_iqtree_starting()

print("EXPERIMENT 6 HIGH SUPPORT")
get_gtm_experiment_6_high_support()
