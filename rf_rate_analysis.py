import sys
import pathlib
import matplotlib.pyplot as plt
import numpy as np
import glob
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
    "100k_seeded": CLUSTER_RUNS + "GTM-experiment6/100k_seeded/METHOD/REPLICATE/",
    # "1000M1": CLUSTER_RUNS + "GTM-experiment6/1000M1/METHOD/REPLICATE/",
    # "aa5K_new": CLUSTER_RUNS + "GTM-experiment6/aa5K_new/METHOD/REPLICATE/",

}

CONSTRAINED_INC_CLUSTER_RUN_MAP = {
    "RandHetCentro10": CLUSTER_RUNS + "PaulRandHetCentro10/Constrained-INC/SUBSET/METHOD/REPLICATE/",
}

GTM_CLUSTER_RUN_MAP = {
    "RandHetCentro10": CLUSTER_RUNS + "PaulRandHetCentro10/GTM/SUBSET/METHOD/REPLICATE/",
}

NJMERGE2_CLUSTER_RUN_MAP = {
    "RandHetCentro10": CLUSTER_RUNS + "PaulRandHetCentro10/NJMerge2/SUBSET/METHOD/REPLICATE/",
}


MODEL_TREE_MAP = {
    "nt10K": PROJECTS_INPUT + "nt10K/huge.REPLICATE.sim.tree",
    "nt78K": PROJECTS_INPUT + "nt78K/huge.REPLICATE.sim.tree",
    "100k_seeded": PROJECTS_INPUT + "100k_seeded/REPLICATE/model/true.tt",
    "1000M1-L": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
    "1000M1-R": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
    "1000M1": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
    "1000S4": PROJECTS_INPUT + "SATe/1000S4/REPLICATE/rose.mt",
    "RNASim1000": PROJECTS_INPUT + "PASTA/RNASim/1000/REPLICATE/model/true.tt",
    "1000M1_HF": PROJECTS_INPUT + "Vlad/UnalignFragTree/high_frag/1000M1/REPLICATE/true_tree.tre",
    "aa5K_new": PROJECTS_INPUT + "aa5K_new/REPLICATE.sim.trim.tree",
    "RandHetCentro10": PROJECTS_INPUT + "Paul/RandHetCentro10/REPLICATE/true-tree.tre",
    "PaulRandHetCentro10": PROJECTS_INPUT + "Paul/RandHetCentro10/REPLICATE/true-tree.tre",
}

REPLICATE_MAP = {
    "nt10K": ["1", "2", "3", "4", "5"],
    "nt78K": ["1", "2", "3", "4", "5"],
    "100k_seeded": ["0", "1", "2", "3", "4"],
    "1000M1-L": ["R0", "R1", "R2", "R3", "R4"],
    "1000M1-R": ["R0", "R1", "R2", "R3", "R4"],
    "1000M1": ["R0", "R1", "R2", "R3", "R4"],
    "1000S4": ["R0", "R1", "R2", "R3", "R4"],
    "RNASim1000": ["1", "2", "3", "4", "5"],
    "1000M1_HF": ["R0", "R1", "R2", "R3", "R4"],
    "aa5K_new": ["COG1028", "COG1309", "COG2814", "COG438", "COG583", "COG596", "COG642"],
    "RandHetCentro10": ["R01", "R02", "R03", "R04", "R05"], #, "R06", "R07", "R08", "R09", "R10"],
    "PaulRandHetCentro10": ["R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09", "R10"],
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
    "iqtree2": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment4/1000M1/250/REPLICATE/output/guide.treefile",
}

BASE_1000M1_ERR_MAP = {
    # "raxmlng": "/projects/tallis/minhyuk2/cluster_runs/individual_test/1000M1RAxML/REPLICATE/raxmlng.raxml.log",
    # "fasttree": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment4/1000M1/250/REPLICATE/errors/fasttree.err",
    # "iqtree2": "/projects/tallis/minhyuk2/cluster_runs/GTM-experiment4/1000M1/250/REPLICATE/output/guide.log",
}

BASE_RANDHETCENTRO10_MAP = {
    "fasttree": "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/fasttree.tre",
    "iqtree2": "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/iqtree+GTR+G/iqtree-result.treefile",
    "raxmlng": "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/raxml+GTR+G/raxml+GTR+G.raxml.bestTree",
}

BASE_aa5K_new_MAP = {
    "raxmlng": "/projects/tallis/minhyuk2/cluster_runs/aa5K_new/RAxML-ng/REPLICATE/raxmlng_output.raxml.lastTree.TMP",
    "iqtree2": "/projects/tallis/minhyuk2/cluster_runs/aa5K_new/IQTree2/REPLICATE/iqtree2_output.treefile",
}

BASE_100k_seeded_MAP = {
    "fasttree": CLUSTER_RUNS + "GTM-experiment6/100k_seeded/IQTree2/REPLICATE/output/fasttree.out",
}

STARTING_1000M1_MAP = {
    "fasttree": "/projects/tallis/minhyuk2/cluster_runs/1000M1_HF/CreateConstraintTrees/FastTree/SUBSET/REPLICATE/output/fasttree.out",
    "clustalomega": "/projects/tallis/minhyuk2/cluster_runs/1000M1_HF/CreateConstraintTrees/ClustalOmega/120/R0/output/clustalomega.tree",
    "mafft": "/projects/tallis/minhyuk2/cluster_runs/1000M1_HF/CreateConstraintTrees/MAFFT/SUBSET/REPLICATE/output/mafft.tree",
    "iqtree2": "/projects/tallis/minhyuk2/cluster_runs/1000M1_HF/IQTree/REPLICATE/iqtree-result.treefile",
    "raxmlng": "/projects/tallis/minhyuk2/cluster_runs/1000M1_HF/RAxML-ng/REPLICATE/raxmlng-result.raxml.lastTree.TMP",
}

STARTING_RANDHETCENTRO10_MAP = {
    "fasttree": "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/fasttree.tre",
    "clustalomega": "/projects/tallis/minhyuk2/cluster_runs/PaulRandHetCentro10/CreateConstraintTrees/ClustalOmega/SUBSET/REPLICATE/output/clustalomega.tree",
    "mafft": "/projects/tallis/minhyuk2/cluster_runs/PaulRandHetCentro10/CreateConstraintTrees/MAFFT/SUBSET/REPLICATE/output/mafft.tree",
    "iqtree2": "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/iqtree+GTR+G/iqtree-result.treefile",
    "raxmlng": "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/raxml+GTR+G/raxml+GTR+G.raxml.bestTree",
}

STARTING_1000S4_MAP = {
    "fasttree": CLUSTER_RUNS + "1000S4/CreateConstraintTrees/SUBSET/FastTree/REPLICATE/output/fasttree.out",
    "iqtree2": CLUSTER_RUNS + "1000S4/CreateConstraintTrees/SUBSET/IQTree2/REPLICATE/output/iqtree-full.treefile",
    "mafft": CLUSTER_RUNS + "1000S4/CreateConstraintTrees/SUBSET/MAFFT/REPLICATE/output/mafft.tree",
}

STARTING_1000M1_HF_MAP = {
    "fasttree": CLUSTER_RUNS + "1000M1_HF/CreateConstraintTrees/FastTree/SUBSET/REPLICATE/output/fasttree.out",
    "iqtree2": CLUSTER_RUNS + "1000M1_HF/IQTree/REPLICATE/iqtree-result.treefile",
}

STARTING_RNASim1000_MAP = {
    "fasttree": CLUSTER_RUNS + "RNASim1000/CreateConstraintTrees/SUBSET/FastTree/REPLICATE/output/fasttree.out",
    "iqtree2": CLUSTER_RUNS + "RNASim1000/CreateConstraintTrees/SUBSET/IQTree2/REPLICATE/output/iqtree-full.treefile",
    "mafft": CLUSTER_RUNS + "RNASim1000/CreateConstraintTrees/SUBSET/MAFFT/REPLICATE/output/mafft.tree",
}

DTM_CONSTRAINED_INC_CLUSTER_RUN_MAP = {
    "fasttree_iqtree": CLUSTER_RUNS + "DATASET/Constrained-INC/SUBSET/fasttree_iqtree/REPLICATE/DIST",
    "fasttree_fasttree": CLUSTER_RUNS + "DATASET/Constrained-INC/SUBSET/fasttree_fasttree/REPLICATE/DIST",
    "mafft": CLUSTER_RUNS + "DATASET/Constrained-INC/SUBSET/MAFFT/REPLICATE/DIST",
    "mafft_mafft": CLUSTER_RUNS + "DATASET/Constrained-INC/SUBSET/mafft_mafft/REPLICATE/DIST",
    "iqtree2": CLUSTER_RUNS + "DATASET/Constrained-INC/SUBSET/iqtree/REPLICATE/DIST",
}

DTM_TREEMERGE_CLUSTER_RUN_MAP = {
    "fasttree_iqtree": CLUSTER_RUNS + "DATASET/TreeMerge/SUBSET/fasttree_iqtree/REPLICATE/DIST",
    "fasttree_fasttree": CLUSTER_RUNS + "DATASET/TreeMerge/SUBSET/fasttree_fasttree/REPLICATE/DIST",
    "mafft": CLUSTER_RUNS + "DATASET/TreeMerge/SUBSET/MAFFT/REPLICATE/DIST",
    "mafft_mafft": CLUSTER_RUNS + "DATASET/TreeMerge/SUBSET/mafft_mafft/REPLICATE/DIST",
    "iqtree2": CLUSTER_RUNS + "DATASET/TreeMerge/SUBSET/iqtree/REPLICATE/DIST",
}

DTM_GTM_CLUSTER_RUN_MAP = {
    "fasttree_iqtree": CLUSTER_RUNS + "DATASET/GTM/SUBSET/fasttree_iqtree/REPLICATE/branch_length.",
    "fasttree_fasttree": CLUSTER_RUNS + "DATASET/GTM/SUBSET/fasttree_fasttree/REPLICATE/branch_length.",
    "mafft": CLUSTER_RUNS + "DATASET/GTM/SUBSET/MAFFT/REPLICATE/branch_length.",
    "mafft_mafft": CLUSTER_RUNS + "DATASET/GTM/SUBSET/mafft_mafft/REPLICATE/branch_length.",
    "iqtree2": CLUSTER_RUNS + "DATASET/GTM/SUBSET/iqtree/REPLICATE/branch_length.",
}

NEW_OLD_TREEMERGE_TIME_MEMORY_CLUSTER_RUN_MAP = {
    "default": CLUSTER_RUNS + "DATASET/NewVSOldTreeMerge/METHOD/SUBSET/DTM/REPLICATE/output/treemerge.tree",
}

NEW_TREEMERGE_TIME_MEMORY_CLUSTER_RUN_MAP = {
    "default": CLUSTER_RUNS + "DATASET/NewTreeMerge/DTM/METHOD/SUBSET/REPLICATE/output/treemerge_DIST.tree",
}

DTM_TIME_MEMORY_CLUSTER_RUN_MAP = {
    "default": CLUSTER_RUNS + "DATASET/DTM/METHOD/SUBSET/REPLICATE/DIST",
}

PAIRWISE_EXPERIMENT_CLUSTER_RUN_MAP = {
    "default": CLUSTER_RUNS + "PairMergeExperiment/METHOD/DATASET/REPLICATE/DIST",
}


def get_time_memory_from_file(log_file):
    current_time = None
    current_memory = None
    with open(log_file, "r") as f:
        for current_line in f:
            if "Elapsed (wall clock) time (h:mm:ss or m:ss):" in current_line:
                current_line_arr = current_line.split(":")
                hours = 0.0
                minutes = 0.0
                seconds = 0.0
                if len(current_line_arr) == 6:
                    # we are m:ss.xx
                    minutes = int(current_line_arr[4])
                    seconds = int(current_line_arr[5].split(".")[0])
                elif len(current_line_arr) == 7:
                    # we are h:mm:ss.xx
                    hours = int(current_line_arr[4])
                    minutes = int(current_line_arr[5])
                    seconds = int(current_line_arr[6].split(".")[0])
                current_time = (hours * 60 * 60) + (minutes * 60) + seconds
            if "Maximum resident set size (kbytes):" in current_line:
                current_line_arr = current_line.split(":")
                current_memory = int(current_line_arr[1])
    return current_time,current_memory


def get_time_memory(cluster_run_path, method, dist):
    current_time = None
    current_memory = None
    log_file = None

    if(dist == "node."):
        log_file = cluster_run_path + "/" + method + "node.err"
    elif(dist == "branch_length."):
        log_file = cluster_run_path + "/" + method + "brlen.err"
    else:
        log_file = cluster_run_path + "/" + method + ".err"
    return get_time_memory_from_file(log_file)



def get_time_starting(cluster_run_path):
    starting_tree_file = cluster_run_path + "errors/fasttree.err"
    return get_time_fasttree(starting_tree_file)

def get_time_fasttree(starting_tree_file):
    with open(starting_tree_file, "r") as f:
        for line in f:
            if "Total time:" in line:
                time_line = float(line.split(" ")[2])
                return time_line # seconds
                # return time_line * (1/60) * (1/60) # hours


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
            fp_rate += (fp / ei2)
            # print(f'fn_rate: {fn / ei1} and fp_rate: {fp / ei2}')
        elif(constraint_method == "RAxML-ng"):
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, cluster_run_path + "output/sequence_partition_" + str(i) + ".raxml.bestTree")
            fn_rate += (fn / ei1)
            fp_rate += (fp / ei2)
            # print(f'fn_rate: {fn / ei1} and fp_rate: {fp / ei2}')
        elif(constraint_method == "IQTree2"):
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, cluster_run_path + "output/sequence_partition_" + str(i) + ".treefile")
            fn_rate += (fn / ei1)
            fp_rate += (fp / ei2)
            # print(f'fn_rate: {fn / ei1} and fp_rate: {fp / ei2}')
        # if(fp != fn):
        #     print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
        #     print("error: fp != fn")

    fn_rate /= num_clusters
    fp_rate /= num_clusters
    return {
        "fn": fn_rate,
        "fp": fp_rate,
    }

def get_paul_constraint_accuracies(starting_tree, subset_size, replicate):
    base_path = "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/".replace("REPLICATE", replicate)
    fn_rate = 0.0
    fp_rate = 0.0

    if(starting_tree == "fasttree"):
        base_path += "fasttree-centro-"
    elif(starting_tree == "iqtree"):
        base_path += "iqtree+GTR+G-centro-"
    base_path += str(subset_size)
    base_path += "/centro-" + str(subset_size) + "-CLUSTERNUM.treefile"
    num_clusters = -42
    for i in range(100):
        current_tree = base_path.replace("CLUSTERNUM", str(i))
        if not Path(current_tree).is_file():
            num_clusters = i - 1
            break
    for i in range(num_clusters):
        current_tree = base_path.replace("CLUSTERNUM", str(i))
        current_model_tree = current_tree.replace(".treefile", "-base.treefile")
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
        fn_rate += (fn / ei1)
        fp_rate += (fp / ei2)
    fn_rate /= num_clusters
    fp_rate /= num_clusters
    return {
        "fn": fn_rate,
        "fp": fp_rate,
    }


def get_time_iqtree(iqtree_path):
    with open(iqtree_path, "r") as f:
        for line in f:
            # if "Total wall-clock time used:" in line:
            if "Total CPU time used:" in line:
                time_line = float(line.split(" ")[4])
                return time_line
                # return time_line * (1/60) * (1/60) # hours

def get_memory_iqtree(iqtree_path):
    with open(iqtree_path, "r") as f:
        for line in f:
            if "required" in line:
                memory_line = float(line.split(" ")[1])
                return memory_line * 1000

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
                cumulative_fp_rate += (fp / ei2)
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
            cumulative_fp_rate += (fp / ei2)
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
            cumulative_fp_rate += (fp / ei2)
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
            cumulative_fasttree_fp_rate += (fp / ei2)
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
            cumulative_internal_fp_rate += (fp / ei2)
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
                cumulative_fp_rate += (fp / ei2)
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
                cumulative_fp_rate += (fp / ei2)
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
            # cumulative_fp_rate += (fp / ei2)

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
            cumulative_fp_rate += (fp / ei2)
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
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method).replace("aa5K_new", "aa5K_new_" + str(support_threshold) + "_decompose")
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
                cumulative_fp_rate += (fp / ei2)
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fn rate: " + str(fn / ei1))
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fp rate: " + str(fn / ei1))

                nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                if(fp != fn):
                    print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    print("error: fp != fn")
                cumulative_guide_fn_rate += (fn / ei1)
                cumulative_guide_fp_rate += (fp / ei2)

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
            cumulative_fp_rate += (fp / ei2)

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
                cumulative_fp_rate += (fp / ei2)
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
                cumulative_fp_rate += (fp / ei2)

                nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                if(fp != fn):
                    print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    print("error: fp != fn")
                cumulative_guide_fn_rate += (fn / ei1)
                cumulative_guide_fp_rate += (fp / ei2)

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


def get_gtm_experiment_6_1000M1_support_decompose(support_threshold):
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
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method).replace("1000M1", "1000M1_" + str(support_threshold) + "_decompose")
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
                cumulative_fp_rate += (fp / ei2)
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fn rate: " + str(fn / ei1))
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fp rate: " + str(fn / ei1))

                nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                if(fp != fn):
                    print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    print("error: fp != fn")
                cumulative_guide_fn_rate += (fn / ei1)
                cumulative_guide_fp_rate += (fp / ei2)

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


def get_gtm_experiment_6_aa5K_centroid():
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
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method)
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
                cumulative_fp_rate += (fp / ei2)
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fn rate: " + str(fn / ei1))
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fp rate: " + str(fn / ei1))

                nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                if(fp != fn):
                    print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    print("error: fp != fn")
                cumulative_guide_fn_rate += (fn / ei1)
                cumulative_guide_fp_rate += (fp / ei2)

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


def get_gtm_experiment_6_aa5K_support_decompose(support_threshold):
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
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method).replace("aa5K_new", "aa5K_new_" + str(support_threshold) + "_decompose")
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
                cumulative_fp_rate += (fp / ei2)
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fn rate: " + str(fn / ei1))
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fp rate: " + str(fn / ei1))

                # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                # if(fp != fn):
                    # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    # print("error: fp != fn")
                # cumulative_guide_fn_rate += (fn / ei1)
                # cumulative_guide_fp_rate += (fp / ei2)

                print("dataset: " + dataset + " method: " + constraint_method  + " replicate: " + replicate)
                sys.stdout.flush()
                # constraint_accuracies = get_constraint_accuracies(current_cluster_run_path, constraint_method)
                # cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                # cumulative_constraint_fp_rate += constraint_accuracies["fp"]

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
                # average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)

                print(dataset + "-" + constraint_method + "-gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + constraint_method + "-gtm fp rate: " + str(average_fp_rate))
                # print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                # print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                print("starting tree time: " + str(np.median(time_starting_arr)))
                print("constraints tree time: " + str(np.median(time_constraints_arr)))
                print("guide tree time: " + str(np.median(time_guide_arr)))
                print("merge tree time: " + str(np.median(time_merge_arr)))

def get_gtm_experiment_6_fasttree(support_threshold):
    for dataset in EXPERIMENT_6_CLUSTER_RUN_MAP:
        cumulative_fasttree_fn_rate = 0.0
        cumulative_fasttree_fp_rate = 0.0
        time_fasttree_arr = []
        fasttree_missing_count = 0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", "FastTree").replace("aa5K_new", "aa5K_new_" + str(support_threshold) + "_decompose")
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
            cumulative_fasttree_fp_rate += (fp / ei2)
            print(dataset + " replicate: " + replicate + "fn rate: " + str(fn / ei1))
            print(dataset + " replicate: " + replicate + "fp rate: " + str(fn / ei1))
            time_fasttree_arr.append(get_time_fasttree(current_cluster_run_path + "errors/fasttree.err"))

        average_fasttree_fn_rate = cumulative_fasttree_fn_rate / (len(REPLICATE_MAP[dataset]) - fasttree_missing_count)
        average_fasttree_fp_rate = cumulative_fasttree_fp_rate / (len(REPLICATE_MAP[dataset]) - fasttree_missing_count)
        print(dataset + "-fasttree fn rate: " + str(average_fasttree_fn_rate))
        print(dataset + "-fasttree fp rate: " + str(average_fasttree_fp_rate))
        print("fasttree tree time: " + str(np.median(time_fasttree_arr)))


def get_gtm_experiment_6_aa5K_1000_support_decompose(support_threshold):
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
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method).replace("aa5K_new", "aa5K_new_1000_" + str(support_threshold))
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
                cumulative_fp_rate += (fp / ei2)
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fn rate: " + str(fn / ei1))
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fp rate: " + str(fn / ei1))

                # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                # if(fp != fn):
                    # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    # print("error: fp != fn")
                # cumulative_guide_fn_rate += (fn / ei1)
                # cumulative_guide_fp_rate += (fp / ei2)

                print("dataset: " + dataset + " method: " + constraint_method  + " replicate: " + replicate)
                sys.stdout.flush()
                # constraint_accuracies = get_constraint_accuracies(current_cluster_run_path, constraint_method)
                # cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                # cumulative_constraint_fp_rate += constraint_accuracies["fp"]

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
                # average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)

                print(dataset + "-" + constraint_method + "-gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + constraint_method + "-gtm fp rate: " + str(average_fp_rate))
                # print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                # print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                print("starting tree time: " + str(np.median(time_starting_arr)))
                print("constraints tree time: " + str(np.median(time_constraints_arr)))
                print("guide tree time: " + str(np.median(time_guide_arr)))
                print("merge tree time: " + str(np.median(time_merge_arr)))


def get_gtm_experiment_6_aa5K_1000_FT_guide_support_decompose(support_threshold):
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
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method).replace("aa5K_new", "aa5K_new_1000_FTGuide_" + str(support_threshold))
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
                cumulative_fp_rate += (fp / ei2)
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fn rate: " + str(fn / ei1))
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fp rate: " + str(fn / ei1))

                # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                # if(fp != fn):
                    # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    # print("error: fp != fn")
                # cumulative_guide_fn_rate += (fn / ei1)
                # cumulative_guide_fp_rate += (fp / ei2)

                print("dataset: " + dataset + " method: " + constraint_method  + " replicate: " + replicate)
                sys.stdout.flush()
                # constraint_accuracies = get_constraint_accuracies(current_cluster_run_path, constraint_method)
                # cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                # cumulative_constraint_fp_rate += constraint_accuracies["fp"]

                time_starting_arr.append(get_time_starting(current_cluster_run_path))
                time_constraints_arr.append(get_time_constraints(current_cluster_run_path, method=constraint_method))
                # if(dataset == "1000M1"):
                    # time_guide_arr.append(get_time_iqtree(current_cluster_run_path + "output/iqtree2-complete.out"))
                # else:
                    # time_guide_arr.append(get_time_guide(current_cluster_run_path))
                time_merge_arr.append(get_time_merge(current_cluster_run_path))
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + constraint_method + " is not done yet.")
                print(str(missing_count) + " numbers of replicates are missing.")
            else:
                average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)

                print(dataset + "-" + constraint_method + "-gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + constraint_method + "-gtm fp rate: " + str(average_fp_rate))
                # print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                # print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                print("starting tree time: " + str(np.median(time_starting_arr)))
                print("constraints tree time: " + str(np.median(time_constraints_arr)))
                # print("guide tree time: " + str(np.median(time_guide_arr)))
                print("merge tree time: " + str(np.median(time_merge_arr)))



def get_aa5K_new_base():
    dataset = "aa5K_new"
    for method in ["raxmlng", "iqtree2"]:
        cumulative_fn_rate = 0.0
        cumulative_fp_rate = 0.0
        cumulative_time = 0.0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_tree = BASE_aa5K_new_MAP[method].replace("REPLICATE", replicate)
            if not Path(current_tree).is_file():
                # missing_count += 1
                print(current_tree + " is missing")
                continue
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
            if(fp != fn):
                print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                print("error: fp != fn")
            cumulative_fn_rate += (fn / ei1)
            cumulative_fp_rate += (fp / ei2)

            print(dataset + "-" + method + "-" + replicate + " fn rate: ", (fn / ei1))
            print(dataset + "-" + method + "-" + replicate + " fp rate: ", (fp / ei2))
            # if(method == "raxmlng"):
            #     cumulative_time += get_time_raxmlng(error_path)
            # elif(method == "iqtree2"):
            #     cumulative_time += get_time_iqtree(error_path)
            # elif(method == "fasttree"):
            #     cumulative_time += get_time_fasttree(error_path)

        print(dataset + "-" + method + " fn rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
        print(dataset + "-" + method + " fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
        # print(dataset + "-" + method + " average time: ", str(cumulative_time / len(REPLICATE_MAP[dataset])))


def get_gtm_experiment_6_100k_FT_1000_FT():
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
                current_cluster_run_path = None
                if(constraint_method == "IQTree2"):
                    current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method)
                elif(constraint_method == "FastTree" or constraint_method == "RAxML-ng"):
                    current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("100k_seeded/METHOD", "backups/100k_seeded/" + constraint_method)
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
                cumulative_fp_rate += (fp / ei2)
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fn rate: " + str(fn / ei1))
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fp rate: " + str(fn / ei1))

                # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                # if(fp != fn):
                    # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    # print("error: fp != fn")
                # cumulative_guide_fn_rate += (fn / ei1)
                # cumulative_guide_fp_rate += (fp / ei2)

                print("dataset: " + dataset + " method: " + constraint_method  + " replicate: " + replicate)
                sys.stdout.flush()
                # constraint_accuracies = get_constraint_accuracies(current_cluster_run_path, constraint_method)
                # cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                # cumulative_constraint_fp_rate += constraint_accuracies["fp"]

                # time_starting_arr.append(get_time_starting(current_cluster_run_path))
                # time_constraints_arr.append(get_time_constraints(current_cluster_run_path, method=constraint_method))
                # if(dataset == "1000M1"):
                    # time_guide_arr.append(get_time_iqtree(current_cluster_run_path + "output/iqtree2-complete.out"))
                # else:
                    # time_guide_arr.append(get_time_guide(current_cluster_run_path))
                time_merge_arr.append(get_time_merge(current_cluster_run_path))
            if(missing_count == len(REPLICATE_MAP[dataset])):
                print(dataset + "-" + constraint_method + " is not done yet.")
                print(str(missing_count) + " numbers of replicates are missing.")
            else:
                average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                print(dataset + "-" + constraint_method + "-gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + constraint_method + "-gtm fp rate: " + str(average_fp_rate))
                # print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                # print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                # print("starting tree time: " + str(np.median(time_starting_arr)))
                # print("constraints tree time: " + str(np.median(time_constraints_arr)))
                # print("guide tree time: " + str(np.median(time_guide_arr)))
                # print("merge tree time: " + str(np.median(time_merge_arr)))

def get_100k_seeded_base():
    dataset = "100k_seeded"
    for method in ["fasttree"]:
        cumulative_fn_rate = 0.0
        cumulative_fp_rate = 0.0
        cumulative_time = 0.0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_tree = BASE_100k_seeded_MAP[method].replace("REPLICATE", replicate)
            if not Path(current_tree).is_file():
                # missing_count += 1
                print(current_tree + " is missing")
                continue
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
            if(fp != fn):
                print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                print("error: fp != fn")
            cumulative_fn_rate += (fn / ei1)
            cumulative_fp_rate += (fp / ei2)

            print(dataset + "-" + method + "-" + replicate + " fn rate: ", (fn / ei1))
            print(dataset + "-" + method + "-" + replicate + " fp rate: ", (fp / ei2))
            # if(method == "raxmlng"):
            #     cumulative_time += get_time_raxmlng(error_path)
            # elif(method == "iqtree2"):
            #     cumulative_time += get_time_iqtree(error_path)
            # elif(method == "fasttree"):
            #     cumulative_time += get_time_fasttree(error_path)

        print(dataset + "-" + method + " fn rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
        print(dataset + "-" + method + " fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
        # print(dataset + "-" + method + " average time: ", str(cumulative_time / len(REPLICATE_MAP[dataset])))


def get_gtm_experiment_6_100k_FT_heuristic_1000_FT():
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
                current_cluster_run_path = None
                current_cluster_run_path = EXPERIMENT_6_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", constraint_method).replace("100k_seeded", "100k_seeded_FT_heuristic_1000_FT_GTM")
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
                cumulative_fp_rate += (fp / ei2)
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fn rate: " + str(fn / ei1))
                print(dataset + "-" + constraint_method + " replicate: " + replicate + "-gtm fp rate: " + str(fn / ei1))

                # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_cluster_run_path + "output/guide.treefile")
                # if(fp != fn):
                    # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    # print("error: fp != fn")
                # cumulative_guide_fn_rate += (fn / ei1)
                # cumulative_guide_fp_rate += (fp / ei2)

                # print("dataset: " + dataset + " method: " + constraint_method  + " replicate: " + replicate)
                # sys.stdout.flush()
                # constraint_accuracies = get_constraint_accuracies(current_cluster_run_path, constraint_method)
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
                average_fn_rate = cumulative_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_fp_rate = cumulative_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                print(dataset + "-" + constraint_method + "-gtm fn rate: " + str(average_fn_rate))
                print(dataset + "-" + constraint_method + "-gtm fp rate: " + str(average_fp_rate))
                # print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                # print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                # print("starting tree time: " + str(np.median(time_starting_arr)))
                # print("constraints tree time: " + str(np.median(time_constraints_arr)))
                # print("guide tree time: " + str(np.median(time_guide_arr)))
                # print("merge tree time: " + str(np.median(time_merge_arr)))


def get_paul_constrained_inc():
    for dataset in CONSTRAINED_INC_CLUSTER_RUN_MAP:
        for starting_tree in ["fasttree_iqtree"]: #["iqtree+GTR+G"]: #["iqtree+GTR+G", "fasttree"]:
            for subset_size in ["120", "500"]:
                sys.stdout.flush()
                cumulative_node_dist_fn_rate = 0.0
                cumulative_node_dist_fp_rate = 0.0
                cumulative_brlen_fn_rate = 0.0
                cumulative_brlen_fp_rate = 0.0

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
                    current_cluster_run_path = None
                    current_cluster_run_path = CONSTRAINED_INC_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", starting_tree).replace("SUBSET", subset_size)
                    current_node_dist_tree = current_cluster_run_path + "node_dist."
                    current_brlen_tree = current_cluster_run_path + "branch_length."

                    if not Path(current_node_dist_tree).is_file():
                        print(current_node_dist_tree + " is missing")
                        continue
                    if not Path(current_brlen_tree).is_file():
                        print(current_brlen_tree + " is missing")
                        continue

                    nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_node_dist_tree)
                    if(fp != fn):
                        print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        print("error: fp != fn")
                    cumulative_node_dist_fn_rate += (fn / ei1)
                    cumulative_node_dist_fp_rate += (fp / ei2)
                    nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_brlen_tree)
                    if(fp != fn):
                        print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        print("error: fp != fn")
                    cumulative_brlen_fn_rate += (fn / ei1)
                    cumulative_brlen_fp_rate += (fp / ei2)
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " node dist Constraint-INC fn rate: " + str(fn / ei1))
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " node dist Constraint-INC fp rate: " + str(fn / ei1))
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " brlen Constraint-INC fn rate: " + str(fn / ei1))
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " brlen Constraint-INC fp rate: " + str(fn / ei1))

                    # guide_tree = None
                    # guide_tree = "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/fasttree.tre".replace("REPLICATE", replicate)

                    # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, guide_tree)
                    # if(fp != fn):
                        # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        # print("error: fp != fn")
                    # cumulative_guide_fn_rate += (fn / ei1)
                    # cumulative_guide_fp_rate += (fp / ei2)

                    # constraint_accuracies = get_paul_constraint_accuracies(starting_tree, subset_size, replicate)
                    # cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                    # cumulative_constraint_fp_rate += constraint_accuracies["fp"]
                    # time_constraints_arr.append(get_paul_time_constraints(starting_tree, subset_size, replicate))

                average_node_dist_fn_rate = cumulative_node_dist_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_node_dist_fp_rate = cumulative_node_dist_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_brlen_fn_rate = cumulative_brlen_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_brlen_fp_rate = cumulative_brlen_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                print(dataset + "-" + starting_tree + "-" + subset_size + "node_dist-Constrained-INC fn rate: " + str(average_node_dist_fn_rate))
                print(dataset + "-" + starting_tree + "-" + subset_size + "node_dist-Constrained-INC fp rate: " + str(average_node_dist_fp_rate))
                print(dataset + "-" + starting_tree + "-" + subset_size + "brlen-Constrained-INC fn rate: " + str(average_brlen_fn_rate))
                print(dataset + "-" + starting_tree + "-" + subset_size + "brlen-Constrained-INC fp rate: " + str(average_brlen_fp_rate))
                # print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                # print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                # print("mean constraints tree time: " + str(np.mean(time_constraints_arr)))

def get_paul_gtm():
    for dataset in GTM_CLUSTER_RUN_MAP:
        for starting_tree in ["fasttree_iqtree"]: #["iqtree+GTR+G"]: #["iqtree+GTR+G", "fasttree"]:
            for subset_size in ["120", "500"]:
                sys.stdout.flush()
                cumulative_brlen_fn_rate = 0.0
                cumulative_brlen_fp_rate = 0.0

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
                    current_cluster_run_path = None
                    current_cluster_run_path = GTM_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", starting_tree).replace("SUBSET", subset_size)
                    current_brlen_tree = current_cluster_run_path + "branch_length."

                    if not Path(current_brlen_tree).is_file():
                        print(current_brlen_tree + " is missing")
                        continue

                    nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_brlen_tree)
                    if(fp != fn):
                        print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        print("error: fp != fn")
                    cumulative_brlen_fn_rate += (fn / ei1)
                    cumulative_brlen_fp_rate += (fp / ei2)
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " brlen GTM fn rate: " + str(fn / ei1))
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " brlen GTM fp rate: " + str(fn / ei1))

                    # guide_tree = None
                    # guide_tree = "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/fasttree.tre".replace("REPLICATE", replicate)

                    # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, guide_tree)
                    # if(fp != fn):
                        # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        # print("error: fp != fn")
                    # cumulative_guide_fn_rate += (fn / ei1)
                    # cumulative_guide_fp_rate += (fp / ei2)

                    # constraint_accuracies = get_paul_constraint_accuracies(starting_tree, subset_size, replicate)
                    # cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                    # cumulative_constraint_fp_rate += constraint_accuracies["fp"]
                    # time_constraints_arr.append(get_paul_time_constraints(starting_tree, subset_size, replicate))

                average_brlen_fn_rate = cumulative_brlen_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_brlen_fp_rate = cumulative_brlen_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                print(dataset + "-" + starting_tree + "-" + subset_size + "brlen-GTM fn rate: " + str(average_brlen_fn_rate))
                print(dataset + "-" + starting_tree + "-" + subset_size + "brlen-GTM fp rate: " + str(average_brlen_fp_rate))
                # print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                # print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                # print("mean constraints tree time: " + str(np.mean(time_constraints_arr)))

def get_paul_njmerge2():
    for dataset in NJMERGE2_CLUSTER_RUN_MAP:
        for starting_tree in ["fasttree_iqtree"]: #["iqtree+GTR+G"]: #["iqtree+GTR+G", "fasttree"]:
            for subset_size in ["120", "500"]:
                sys.stdout.flush()
                cumulative_node_dist_fn_rate = 0.0
                cumulative_node_dist_fp_rate = 0.0
                cumulative_brlen_fn_rate = 0.0
                cumulative_brlen_fp_rate = 0.0

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
                    current_cluster_run_path = None
                    current_cluster_run_path = NJMERGE2_CLUSTER_RUN_MAP[dataset].replace("REPLICATE", replicate).replace("METHOD", starting_tree).replace("SUBSET", subset_size)
                    current_node_dist_tree = current_cluster_run_path + "node_dist."
                    current_brlen_tree = current_cluster_run_path + "branch_length."

                    if not Path(current_node_dist_tree).is_file():
                        print(current_node_dist_tree + " is missing")
                        continue
                    if not Path(current_brlen_tree).is_file():
                        print(current_brlen_tree + " is missing")
                        continue

                    nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_node_dist_tree)
                    if(fp != fn):
                        print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        print("error: fp != fn")
                    cumulative_node_dist_fn_rate += (fn / ei1)
                    cumulative_node_dist_fp_rate += (fp / ei2)
                    nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_brlen_tree)
                    if(fp != fn):
                        print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        print("error: fp != fn")
                    cumulative_brlen_fn_rate += (fn / ei1)
                    cumulative_brlen_fp_rate += (fp / ei2)
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " node dist GTM fn rate: " + str(fn / ei1))
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " node dist GTM fp rate: " + str(fn / ei1))
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " brlen GTM fn rate: " + str(fn / ei1))
                    print(dataset + "-" + starting_tree + "-" + subset_size + " replicate: " + replicate + " brlen GTM fp rate: " + str(fn / ei1))

                    # guide_tree = None
                    # guide_tree = "/projects/tallis/minhyuk2/paul_rhc10/random-het-centro10/REPLICATE/fasttree.tre".replace("REPLICATE", replicate)

                    # nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, guide_tree)
                    # if(fp != fn):
                        # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        # print("error: fp != fn")
                    # cumulative_guide_fn_rate += (fn / ei1)
                    # cumulative_guide_fp_rate += (fp / ei2)

                    # constraint_accuracies = get_paul_constraint_accuracies(starting_tree, subset_size, replicate)
                    # cumulative_constraint_fn_rate += constraint_accuracies["fn"]
                    # cumulative_constraint_fp_rate += constraint_accuracies["fp"]
                    # time_constraints_arr.append(get_paul_time_constraints(starting_tree, subset_size, replicate))

                average_node_dist_fn_rate = cumulative_node_dist_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_node_dist_fp_rate = cumulative_node_dist_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_brlen_fn_rate = cumulative_brlen_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                average_brlen_fp_rate = cumulative_brlen_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fn_rate = cumulative_guide_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_guide_fp_rate = cumulative_guide_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fn_rate = cumulative_constraint_fn_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                # average_constraint_fp_rate = cumulative_constraint_fp_rate / (len(REPLICATE_MAP[dataset]) - missing_count)
                print(dataset + "-" + starting_tree + "-" + subset_size + "node_dist-gtm fn rate: " + str(average_node_dist_fn_rate))
                print(dataset + "-" + starting_tree + "-" + subset_size + "node_dist-gtm fp rate: " + str(average_node_dist_fp_rate))
                print(dataset + "-" + starting_tree + "-" + subset_size + "brlen-gtm fn rate: " + str(average_brlen_fn_rate))
                print(dataset + "-" + starting_tree + "-" + subset_size + "brlen-gtm fp rate: " + str(average_brlen_fp_rate))
                # print(dataset + "-guide fn rate: " + str(average_guide_fn_rate))
                # print(dataset + "-guide fp rate: " + str(average_guide_fp_rate))
                # print(dataset + "-constraint fn rate: " + str(average_constraint_fn_rate))
                # print(dataset + "-constraint fp rate: " + str(average_constraint_fp_rate))
                # print("mean constraints tree time: " + str(np.mean(time_constraints_arr)))

def get_randhetcentro10_base():
    dataset = "RandHetCentro10"
    for method in ["fasttree", "iqtree2", "raxmlng"]:
        cumulative_fn_rate = 0.0
        cumulative_fp_rate = 0.0
        cumulative_time = 0.0
        for replicate in REPLICATE_MAP[dataset]:
            current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
            current_tree = BASE_RANDHETCENTRO10_MAP[method].replace("REPLICATE", replicate)
            if not Path(current_tree).is_file():
                # missing_count += 1
                print(current_tree + " is missing")
                continue
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
            if(fp != fn):
                print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                print("error: fp != fn")
            cumulative_fn_rate += (fn / ei1)
            cumulative_fp_rate += (fp / ei2)

            print(dataset + "-" + method + "-" + replicate + " fn rate: ", (fn / ei1))
            print(dataset + "-" + method + "-" + replicate + " fp rate: ", (fp / ei2))
            # if(method == "raxmlng"):
            #     cumulative_time += get_time_raxmlng(error_path)
            # elif(method == "iqtree2"):
            #     cumulative_time += get_time_iqtree(error_path)
            # elif(method == "fasttree"):
            #     cumulative_time += get_time_fasttree(error_path)

        print(dataset + "-" + method + " fn rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
        print(dataset + "-" + method + " fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
        # print(dataset + "-" + method + " average time: ", str(cumulative_time / len(REPLICATE_MAP[dataset])))


def get_randhetcentro10_1000M1_starting():
    dataset_arr = ["RandHetCentro10", "1000M1"]
    for dataset in dataset_arr:
        for method in ["fasttree", "iqtree2", "clustalomega", "mafft", "raxmlng"]:
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_tree = None
                if dataset == "RandHetCentro10":
                    current_tree = STARTING_RANDHETCENTRO10_MAP[method].replace("REPLICATE", replicate).replace("SUBSET", "120")
                elif dataset == "1000M1":
                    current_tree = STARTING_1000M1_MAP[method].replace("REPLICATE", replicate).replace("SUBSET", "120")
                if not Path(current_tree).is_file():
                    # missing_count += 1
                    print(current_tree + " is missing")
                    continue
                nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                # if(fp != fn):
                    # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                    # print("error: fp != fn")
                cumulative_fn_rate += (fn / ei1)
                cumulative_fp_rate += (fp / ei2)

                # print("starting " + dataset + "-" + method + "-" + replicate + " fn rate: ", (fn / ei1))
                # print("starting " + dataset + "-" + method + "-" + replicate + " fp rate: ", (fp / ei2))

            # print("starting " + dataset + "-" + method + " fn rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
            # print("starting " + dataset + "-" + method + " fp rate: ", str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
            print("starting " + dataset + "-" + method + " fn/fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])) + "/" + str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))


def get_DTM_constrained_inc():
    meta_arr = []
    dataset_arr = ["1000M1_HF", "RNASim1000"] # ["PaulRandHetCentro10", "1000M1_HF"]
    for dataset in dataset_arr:
        data_arr = []
        labels = []
        #for method in ["fasttree_iqtree", "IQTree2", "MAFFT", "fasttree_fasttree", "mafft_mafft"]:
        for method in ["fasttree_iqtree"]:
            for subset in ["500"]:
                for dist in ["node."]: #, "branch_length."]:
                    cumulative_fn_rate = 0.0
                    cumulative_fp_rate = 0.0
                    current_data_arr = []
                    time_constraints_arr = []
                    for replicate in REPLICATE_MAP[dataset]:
                        current_time = 0.0
                        current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                        current_tree = DTM_CONSTRAINED_INC_CLUSTER_RUN_MAP[method.lower()].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("fasttree_iqtree", subset).replace("SUBSET", method).replace("DIST", dist)
                        if not Path(current_tree).is_file():
                            # missing_count += 1
                            print(current_tree + " is missing")
                            continue
                        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                        # if(fp != fn):
                            # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                            # print("error: fp != fn")
                        cumulative_fn_rate += (fn / ei1)
                        cumulative_fp_rate += (fp / ei2)

                        # print("Constrained-inc " + dataset + "-" + subset + "-" + method + "-" + dist + "-" + replicate + " fn rate: ", (fn / ei1))
                        # print("Constrained-inc " + dataset + "-" + subset + "-" + method + "-" + dist + "-" + replicate + " fp rate: ", (fp / ei2))
                        current_data_arr.append(rf)
                        current_cluster_run_path = None
                        if(dataset == "1000M1_HF"):
                            current_cluster_run_path = CLUSTER_RUNS + dataset + "/CreateConstraintTrees/FastTree/500/" + replicate + "/"
                            current_time += get_time_fasttree(CLUSTER_RUNS + dataset + "/CreateConstraintTrees/FastTree/500/" + replicate + "/errors/fasttree.err")
                            # print("starting: " + str(current_time))
                            current_time += (get_time_constraints(current_cluster_run_path, method="IQTree2"))
                            # print("constraints: " + str(current_time))
                            current_time += get_time_memory_from_file(CLUSTER_RUNS + dataset + "/Constrained-INC/fasttree_iqtree/500/" + replicate + "/cinc_node.err")[0]
                            # print("cinc: " + str(current_time))
                            time_constraints_arr.append(current_time)
                        elif(dataset == "RNASim1000"):
                            current_cluster_run_path = CLUSTER_RUNS + dataset + "/CreateConstraintTrees/500/FastTree/" + replicate + "/"
                            current_time += get_time_fasttree(CLUSTER_RUNS + dataset + "/CreateConstraintTrees/500/FastTree/" + replicate + "/errors/fasttree.err")
                            # print("starting: " + str(current_time))
                            current_time += (get_time_constraints(current_cluster_run_path, method="IQTree2"))
                            # print("constraints: " + str(current_time))
                            current_time += get_time_memory_from_file(CLUSTER_RUNS + dataset + "/Constrained-INC/fasttree_iqtree/500/" + replicate + "/cinc_node.err")[0]
                            # print("cinc: " + str(current_time))
                            time_constraints_arr.append(current_time)

                    # print("Constrained-inc" + dataset + "-" + subset + "-" + method + "-" + dist + " fn rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
                    # print("Constrained-inc" + dataset + "-" + subset + "-" + method + "-" + dist + " fp rate: ", str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                    print("Constrained-INC " + dataset + "-" + subset + "-" + method + "-" + dist + " fn/fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])) + "/" + str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                    print("Constrained-INC " + dataset + "-" + subset + "-" + method + "-" + dist + " starting + constraints + cinc time: ", np.mean(time_constraints_arr))
                    data_arr.append(current_data_arr)
                    labels.append(dataset + " Constrained-INC " + " " + subset + " " + method + " " + dist)
        meta_arr.append((data_arr, labels))
    return meta_arr

def get_DTM_treemerge():
    meta_arr = []
    dataset_arr = ["PaulRandHetCentro10", "1000M1_HF"]
    for dataset in dataset_arr:
        data_arr = []
        labels = []
        for method in ["fasttree_iqtree", "IQTree2", "MAFFT", "fasttree_fasttree", "mafft_mafft"]:
            for subset in ["120", "500"]:
                for dist in ["node_dist.", "branch_length."]:
                    cumulative_fn_rate = 0.0
                    cumulative_fp_rate = 0.0
                    current_data_arr = []
                    for replicate in REPLICATE_MAP[dataset]:
                        current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                        current_tree = DTM_TREEMERGE_CLUSTER_RUN_MAP[method.lower()].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("METHOD", method).replace("SUBSET", subset).replace("DIST", dist)
                        if not Path(current_tree).is_file():
                            # missing_count += 1
                            print(current_tree + " is missing")
                            continue
                        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                        # if(fp != fn):
                            # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                            # print("error: fp != fn")
                        cumulative_fn_rate += (fn / ei1)
                        cumulative_fp_rate += (fp / ei2)

                        # print("treemerge " + dataset + "-" + subset + "-" + method + "-" + dist + "-" + replicate + " fn rate: ", (fn / ei1))
                        # print("treemerge " + dataset + "-" + subset + "-" + method + "-" + dist + "-" + replicate + " fp rate: ", (fp / ei2))
                        current_data_arr.append(rf)

                    # print("treemerge " + dataset + "-" + subset + "-" + method + "-" + dist + " fn rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
                    # print("treemerge " + dataset + "-" + subset + "-" + method + "-" + dist + " fp rate: ", str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                    print("treemerge " + dataset + "-" + subset + "-" + method + "-" + dist + " fn/fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])) + "/" + str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                    data_arr.append(current_data_arr)
                    labels.append(dataset + " TreeMerge " + " " + subset + " " + method + " " + dist)
        meta_arr.append((data_arr, labels))
    return meta_arr

def get_DTM_GTM():
    meta_arr = []
    GTM_1000M1_HF_time_arr = [0.5047855377197266, 0.4293036460876465, 0.37412357330322266, 0.38379764556884766, 0.3979816436767578]
    GTM_RNASim1000_time_arr = [0.32965874671936035, 0.291461706161499, 0.2992827892303467, 0.2864358425140381, 0.2990751266479492]

    dataset_arr = ["RNASim1000", "1000M1_HF"]
    for dataset in dataset_arr:
        data_arr = []
        labels = []
        for method in ["fasttree_iqtree"]: #, "IQTree2", "MAFFT", "fasttree_fasttree", "mafft_mafft"]:
            for subset in ["500"]:
                cumulative_fn_rate = 0.0
                cumulative_fp_rate = 0.0
                current_data_arr = []
                time_constraints_arr = []
                for rep_index,replicate in enumerate(REPLICATE_MAP[dataset]):
                    current_time = 0.0
                    current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                    current_tree = DTM_GTM_CLUSTER_RUN_MAP[method.lower()].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("METHOD", method).replace("SUBSET", subset)
                    if(dataset == "RNASim1000"):
                        current_tree = DTM_GTM_CLUSTER_RUN_MAP[method.lower()].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("fasttree_iqtree", subset).replace("SUBSET", method).replace("DIST", "branch_length.")
                    if not Path(current_tree).is_file():
                        # missing_count += 1
                        print(current_tree + " is missing")
                        continue
                    nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                    # if(fp != fn):
                        # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        # print("error: fp != fn")
                    cumulative_fn_rate += (fn / ei1)
                    cumulative_fp_rate += (fp / ei2)

                    # print("GTM " + dataset + "-" + subset + "-" + method + "-" + replicate + " fn rate: ", (fn / ei1))
                    # print("GTM " + dataset + "-" + subset + "-" + method + "-" + replicate + " fp rate: ", (fp / ei2))
                    current_data_arr.append(rf)
                    if(dataset == "1000M1_HF"):
                        current_cluster_run_path = CLUSTER_RUNS + dataset + "/CreateConstraintTrees/FastTree/500/" + replicate + "/"
                        current_time += get_time_fasttree(CLUSTER_RUNS + dataset + "/CreateConstraintTrees/FastTree/500/" + replicate + "/errors/fasttree.err")
                        # print("starting: " + str(current_time))
                        current_time += (get_time_constraints(current_cluster_run_path, method="IQTree2"))
                        # print("constraints: " + str(current_time))
                        current_time += GTM_1000M1_HF_time_arr[rep_index]
                        # print("cinc: " + str(current_time))
                        time_constraints_arr.append(current_time)
                    elif(dataset == "RNASim1000"):
                        current_cluster_run_path = CLUSTER_RUNS + dataset + "/CreateConstraintTrees/500/FastTree/" + replicate + "/"
                        current_time += get_time_fasttree(CLUSTER_RUNS + dataset + "/CreateConstraintTrees/500/FastTree/" + replicate + "/errors/fasttree.err")
                        # print("starting: " + str(current_time))
                        current_time += (get_time_constraints(current_cluster_run_path, method="IQTree2"))
                        # print("constraints: " + str(current_time))
                        current_time += GTM_RNASim1000_time_arr[rep_index]
                        # print("cinc: " + str(current_time))
                        time_constraints_arr.append(current_time)

                # print("GTM " + dataset + "-" + subset + "-" + method + " fn rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
                # print("GTM " + dataset + "-" + subset + "-" + method + " fp rate: ", str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                print("GTM " + dataset + "-" + subset + "-" + method + " fn/fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])) + "/" + str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                print("GTM " + dataset + "-" + subset + "-" + method + " total time: ", np.mean(time_constraints_arr))
                data_arr.append(current_data_arr)
                labels.append(dataset + " GTM " + " " + subset + " " + method)
        meta_arr.append((data_arr, labels))
    return meta_arr


def get_pairmerge_500_results():
    meta_arr = []
    dataset_arr = ["1000M1_HF", "RNASim1000"] # ["1000S4", "1000M1_HF", "RandHetCentro10", "RNASim1000"]
    for dataset in dataset_arr:
        data_arr = []
        labels = []
        for method in ["GTM", "Constrained-INC", "NJMerge2"]:
            for dist in ["induced_FT_node", "IQTree_node"]:
                cumulative_fn_rate = 0.0
                cumulative_fp_rate = 0.0
                current_data_arr = []
                for replicate in REPLICATE_MAP[dataset]:
                    current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                    current_tree = PAIRWISE_EXPERIMENT_CLUSTER_RUN_MAP["default"].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("METHOD", method).replace("DIST", dist)
                    if not Path(current_tree).is_file():
                        # missing_count += 1
                        print(current_tree + " is missing")
                        continue
                    nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                    # if(fp != fn):
                        # print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
                        # print("error: fp != fn")
                    cumulative_fn_rate += (fn / ei1)
                    cumulative_fp_rate += (fp / ei2)

                    # print("treemerge " + dataset + "-" + subset + "-" + method + "-" + dist + "-" + replicate + " fn rate: ", (fn / ei1))
                    # print("treemerge " + dataset + "-" + subset + "-" + method + "-" + dist + "-" + replicate + " fp rate: ", (fp / ei2))
                    print(dataset + "-" + method + "-" + dist + "-" + replicate + " fn rate: ", str(fn / ei1))
                    print(dataset + "-" + method + "-" + dist + "-" + replicate + " fp rate: ", str(fp / ei2))
                    current_data_arr.append(rf)

                print(dataset + "-" + method + "-" + dist + " fn rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])))
                print(dataset + "-" + method + "-" + dist + " fp rate: ", str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                data_arr.append(current_data_arr)
                if(dist == "induced_FT_node"):
                    labels.append(dataset + " " + method + " " + "induced FT")
                elif(dist == "IQTree_node"):
                    labels.append(dataset + " " + method + " " + "IQTree")

        meta_arr.append((data_arr, labels))
    return meta_arr

def get_RNASim1000_1000S4_starting():
    dataset_arr = ["RNASim1000", "1000S4"]
    for dataset in dataset_arr:
        for method in ["fasttree", "iqtree2", "mafft"]:
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0
            cumulative_time = 0.0
            cumulative_memory = 0.0
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_tree = None
                if dataset == "RNASim1000":
                    current_tree = STARTING_RNASim1000_MAP[method].replace("REPLICATE", replicate).replace("SUBSET", "120")
                elif dataset == "1000S4":
                    current_tree = STARTING_1000S4_MAP[method].replace("REPLICATE", replicate).replace("SUBSET", "120")
                if not Path(current_tree).is_file():
                    # missing_count += 1
                    print(current_tree + " is missing")
                    continue
                nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                cumulative_fn_rate += (fn / ei1)
                cumulative_fp_rate += (fp / ei2)

                current_time = None
                current_memory = None
                if(method == "fasttree"):
                    current_time,current_memory = get_time_memory_from_file(current_tree.split("output")[0] + "errors/fasttree.err")
                elif(method == "iqtree2"):
                    current_time,current_memory = get_time_memory_from_file(current_tree.split("output")[0] + "errors/iqtree2-complete.err")
                if(method == "mafft"):
                    current_time,current_memory = get_time_memory_from_file(current_tree.split("output")[0] + "errors/mafft.err")
                cumulative_time += current_time
                cumulative_memory += current_memory

            print("starting " + dataset + "-" + method + " fn/fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])) + "/" + str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
            print("starting " + dataset + "-" + method + " time: ", str(cumulative_time / len(REPLICATE_MAP[dataset])))
            print("starting " + dataset + "-" + method + " memory: ", str(cumulative_memory /  len(REPLICATE_MAP[dataset])))

def get_DTM_time_memory(DTM):
    meta_arr = []
    dataset_arr = ["RNASim1000", "1000M1_HF"] #["RNASim1000", "1000S4"]
    for dataset in dataset_arr:
        data_arr = []
        labels = []
        for method in ["IQTree2", "fasttree_iqtree", "MAFFT", "fasttree_fasttree", "mafft_mafft"]:
            for subset in ["120", "500"]:
                for dist in ["node.", "branch_length."]:
                    cumulative_fn_rate = 0.0
                    cumulative_fp_rate = 0.0
                    cumulative_time = 0.0
                    cumulative_memory = 0.0
                    current_data_arr = []
                    for replicate in REPLICATE_MAP[dataset]:
                        current_max_memory = 0.0
                        current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                        if(DTM == "NewTreeMerge"):
                            if(dist == "node."):
                                current_tree = DTM_TIME_MEMORY_CLUSTER_RUN_MAP[DTM].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("METHOD", method).replace("SUBSET", subset).replace("DIST", "node").replace("DTM", DTM)
                            elif(dist == "branch_length."):
                                current_tree = DTM_TIME_MEMORY_CLUSTER_RUN_MAP[DTM].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("METHOD", method).replace("SUBSET", subset).replace("DIST", "brlen").replace("DTM", DTM)
                        else:
                            current_tree = DTM_TIME_MEMORY_CLUSTER_RUN_MAP["default"].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("METHOD", method).replace("SUBSET", subset).replace("DIST", dist).replace("DTM", DTM)
                        if not Path(current_tree).is_file():
                            # missing_count += 1
                            print(current_tree + " is missing")
                            continue
                        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                        cumulative_fn_rate += (fn / ei1)
                        cumulative_fp_rate += (fp / ei2)
                        current_data_arr.append(rf)

                        method_name_err = None
                        if(DTM == "GTM"):
                            method_name_err = "gtm_"
                        elif(DTM == "Constrained-INC"):
                            method_name_err = "cinc_"
                        elif("TreeMerge" in DTM):
                            method_name_err = "treemerge_"

                        if("TreeMerge" in DTM):
                            current_time,current_memory = get_time_memory_from_file(current_tree.replace("output", "errors").replace(".tree", ".err"))
                            cumulative_time += current_time
                            current_max_memory = max(current_max_memory, current_memory)
                            current_time,current_memory = get_time_memory_from_file(current_tree.replace("output", "errors").split("treemerge_")[0] + "merge.err")
                            cumulative_time += current_time
                            current_max_memory = max(current_max_memory, current_memory)
                            current_time,current_memory = get_time_memory_from_file(current_tree.replace("output", "errors").split("treemerge_")[0] + "setup_merger.err")
                            cumulative_time += current_time
                            current_max_memory = max(current_max_memory, current_memory)
                            current_time,current_memory = get_time_memory_from_file(current_tree.replace("output", "errors").split("treemerge_")[0] + "run_merger.err")
                            cumulative_time += current_time
                            current_max_memory = max(current_max_memory, current_memory)
                        else:
                            current_time,current_memory = get_time_memory(current_tree[:-len(dist)], method_name_err, dist)
                            cumulative_time += current_time
                            current_max_memory = max(current_max_memory, current_memory)
                        cumulative_memory += current_max_memory

                    print(DTM + " " + dataset + "-" + subset + "-" + method + "-" + dist + " fn/fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])) + "/" + str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                    # print(DTM + " " + dataset + "-" + subset + "-" + method + "-" + dist + " time: ", str(cumulative_time / len(REPLICATE_MAP[dataset])))
                    # print(DTM + " " + dataset + "-" + subset + "-" + method + "-" + dist + " memory: ", str(cumulative_memory / len(REPLICATE_MAP[dataset])))
                    data_arr.append(current_data_arr)
                    labels.append(dataset + " DTM " + " " + subset + " " + method + " " + dist)
        meta_arr.append((data_arr, labels))
    return meta_arr


def get_new_vs_old_treemerge(DTM):
    meta_arr = []
    dataset_arr = ["RNASim1000", "1000M1_HF"]
    for dataset in dataset_arr:
        data_arr = []
        labels = []
        for method in ["fasttree_fasttree"]:
            for subset in ["500"]:
                for dist in ["notttt itt"]:
                    cumulative_fn_rate = 0.0
                    cumulative_fp_rate = 0.0
                    cumulative_memory = 0.0
                    current_data_arr = []
                    time_arr = []
                    for replicate in REPLICATE_MAP[dataset]:
                        cumulative_time = 0.0
                        current_max_memory = 0.0
                        current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                        current_tree = NEW_OLD_TREEMERGE_TIME_MEMORY_CLUSTER_RUN_MAP["default"].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("METHOD", method).replace("SUBSET", subset).replace("DIST", dist).replace("DTM", DTM)
                        if not Path(current_tree).is_file():
                            # missing_count += 1
                            print(current_tree + " is missing")
                            continue
                        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                        cumulative_fn_rate += (fn / ei1)
                        cumulative_fp_rate += (fp / ei2)
                        current_data_arr.append(rf)

                        method_name_err = None
                        method_name_err = "treemerge"
                        current_time,current_memory = get_time_memory(current_tree[:-len("output/treemerge.tree")] + "errors/", method_name_err, dist)
                        cumulative_time += current_time
                        current_max_memory = max(current_max_memory, current_memory)
                        if(DTM != "OldTreeMerge"):
                            current_time,current_memory = get_time_memory(current_tree[:-len("output/treemerge.tree")] + "errors/", "merge", dist)
                            cumulative_time += current_time
                            current_max_memory = max(current_max_memory, current_memory)
                            current_time,current_memory = get_time_memory(current_tree[:-len("output/treemerge.tree")] + "errors/", "setup_merger", dist)
                            cumulative_time += current_time
                            current_max_memory = max(current_max_memory, current_memory)
                            current_time,current_memory = get_time_memory(current_tree[:-len("output/treemerge.tree")] + "errors/", "run_merger", dist)
                            cumulative_time += current_time
                            current_max_memory = max(current_max_memory, current_memory)
                        cumulative_memory += current_max_memory
                        time_arr.append(cumulative_time)

                    print(time_arr)
                    print(DTM + " " + dataset + "-" + subset + "-" + method + "-" + dist + " fn/fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])) + "/" + str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                    print(DTM + " " + dataset + "-" + subset + "-" + method + "-" + dist + " time: ", np.mean(time_arr))
                    print(DTM + " " + dataset + "-" + subset + "-" + method + "-" + dist + " memory: ", str(cumulative_memory / len(REPLICATE_MAP[dataset])))
                    data_arr.append(current_data_arr)
                    labels.append(dataset + " DTM " + " " + subset + " " + method + " " + dist)
        meta_arr.append((data_arr, labels))
    return meta_arr


def create_box_plot(data_arr, raw_labels, fig_name, colors):
    labels = [" ".join(element.split()[1:]) for element in raw_labels]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    medianprops = {
        "linewidth": 2.0,
        "color": "red",
    }
    boxplot = ax.boxplot(data_arr, patch_artist=True, medianprops=medianprops)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("RF rate")
    ax.set_ylim((0.0, 0.60))
    for patches,color in zip(boxplot["boxes"], colors):
        patches.set_facecolor(color)

    # ax.get_xaxis().tick_bottom()
    # ax.get_yaxis().tick_left()
    fig.savefig(fig_name, bbox_inches="tight")

def pair_merge_create_box_plots(meta_arr):
    for data_arr,labels in meta_arr:
        create_box_plot(data_arr, labels, "./figures/pair_merge_500/" + labels[0].split()[0], ["gray", "white", "gray", "white", "gray", "white"])

def DTM_treemerge_create_box_plots(meta_arr):
    for data_arr,labels in meta_arr:
        create_box_plot(data_arr, labels, "./figures/iqtree_guide/treemerge/" + labels[0].split()[0], ["gray", "white"] * 6)

def DTM_constrained_inc_create_box_plots(meta_arr):
    for data_arr,labels in meta_arr:
        create_box_plot(data_arr, labels, "./figures/iqtree_guide/constrained_inc/" + labels[0].split()[0], ["gray", "white"] * 6)

def DTM_GTM_create_box_plots(meta_arr):
    for data_arr,labels in meta_arr:
        create_box_plot(data_arr, labels, "./figures/iqtree_guide/GTM/" + labels[0].split()[0], ["gray", "white"] * 4)


def get_new_treemerge_time_memory(DTM):
    meta_arr = []
    dataset_arr = ["1000M1_HF", "RNASim1000"]
    for dataset in dataset_arr:
        data_arr = []
        labels = []
        for method in ["IQTree2", "fasttree_fasttree", "mafft_mafft"]: #, "fasttree_fasttree", "mafft_mafft"]:
            for subset in ["500"]: # ["120", "500"]:
                for dist in ["node"]: # ["node", "brlen"]:
                    cumulative_fn_rate = 0.0
                    cumulative_fp_rate = 0.0
                    cumulative_memory = 0.0
                    current_data_arr = []
                    time_arr = []
                    memory_arr = []
                    for replicate in REPLICATE_MAP[dataset]:
                        cumulative_time = 0.0
                        current_max_memory = 0.0
                        current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                        current_tree = None
                        current_tree = NEW_TREEMERGE_TIME_MEMORY_CLUSTER_RUN_MAP["default"].replace("DATASET", dataset).replace("REPLICATE", replicate).replace("METHOD", method).replace("SUBSET", subset).replace("DIST", dist).replace("DTM", DTM)
                        if not Path(current_tree).is_file():
                            # missing_count += 1
                            print(current_tree + " is missing")
                            continue
                        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                        cumulative_fn_rate += (fn / ei1)
                        cumulative_fp_rate += (fp / ei2)
                        current_data_arr.append(rf)
                        folder_method_name = ""
                        err_file_name = ""
                        if(method == "IQTree2"):
                            folder_method_name = "IQTree2"
                        elif(method == "fasttree_fasttree"):
                            folder_method_name = "FastTree"
                        elif(method == "mafft_mafft"):
                            folder_method_name = "MAFFT"

                        if(dataset == "1000M1_HF"):
                            current_cluster_run_path = CLUSTER_RUNS + dataset + "/CreateConstraintTrees/" + folder_method_name + "/500/" + replicate + "/"
                            if(method == "fasttree_fasttree"):
                                cumulative_time += get_time_fasttree(CLUSTER_RUNS + dataset + "/CreateConstraintTrees/" + folder_method_name +  "/500/" + replicate + "/errors/fasttree.err")
                            elif(method == "IQTree2"):
                                cumulative_time += get_time_iqtree(CLUSTER_RUNS + dataset + "/IQTree/" + replicate + "/iqtree-result.log")
                            # elif(method == "mafft_mafft"):
                                # current_time += get_time_iqtree(CLUSTER_RUNS + dataset + "/MAFFT/" + replicate + "/iqtree-result.log")
                            # print("starting: " + str(current_time))
                            cumulative_time += (get_time_constraints(current_cluster_run_path, method="IQTree2"))
                            # print("constraints: " + str(current_time))
                        elif(dataset == "RNASim1000"):
                            current_cluster_run_path = CLUSTER_RUNS + dataset + "/CreateConstraintTrees/500/" + folder_method_name + "/" + replicate + "/"
                            if(method == "fasttree_fasttree"):
                                cumulative_time += get_time_fasttree(CLUSTER_RUNS + dataset + "/CreateConstraintTrees/500/" + folder_method_name + "/" + replicate + "/errors/fasttree.err")
                            elif(method == "IQTree2"):
                                cumulative_time += get_time_iqtree(CLUSTER_RUNS + dataset + "/CreateConstraintTrees/500/" + folder_method_name + "/" + replicate + "/output/iqtree-full.log")
                            # elif(method == "mafft_mafft"):
                                # current_time += get_time_iqtree(CLUSTER_RUNS + dataset + "/CreateConstraintTrees/500/" + folder_method_name + "/" + replicate + "/output/iqtree-full.log")
                            # print("starting: " + str(current_time))
                            # print("starting: " + str(current_time))
                            cumulative_time += (get_time_constraints(current_cluster_run_path, method="IQTree2"))


                        method_name_err = "treemerge_"

                        current_time,current_memory = get_time_memory_from_file(current_tree.replace("output", "errors").replace(".tree", ".err"))
                        # print("treemerge time: " + str(current_time))
                        cumulative_time += current_time
                        # print("CURRENT TIME: " + str(cumulative_time))
                        current_max_memory = max(current_max_memory, current_memory)
                        current_time,current_memory = get_time_memory_from_file(current_tree.replace("output", "errors").split("treemerge_")[0] + "merge.err")
                        # print("merge time: " + str(current_time))
                        cumulative_time += current_time
                        # print("CURRENT TIME: " + str(cumulative_time))
                        current_max_memory = max(current_max_memory, current_memory)
                        # current_time,current_memory = get_time_memory_from_file(current_tree.replace("output", "errors").split("treemerge_")[0] + "setup_merger.err")
                        # cumulative_time += current_time
                        current_max_memory = max(current_max_memory, current_memory)
                        current_time,current_memory = get_time_memory_from_file(current_tree.replace("output", "errors").split("treemerge_")[0] + "run_merger.err")
                        # print("run_merger time: " + str(current_time))
                        cumulative_time += current_time
                        # print("CURRENT TIME: " + str(cumulative_time))
                        current_max_memory = max(current_max_memory, current_memory)
                        iqtree_time = 0.0

                        for iqtree_file in glob.glob(current_tree.split("output")[0] + "/output/njmergepair-sequence_partition_*-and-sequence_partition_*-starting.tree.log"):
                            iqtree_time = max(iqtree_time, get_time_iqtree(iqtree_file))
                            # print("current_iqtree_time : ", str(get_time_iqtree(iqtree_file)))
                            current_max_memory = max(current_max_memory, get_memory_iqtree(iqtree_file))
                        # print("iqtree time: " + str(iqtree_time))
                        cumulative_time += iqtree_time

                        # print("CURRENT TIME: " + str(cumulative_time))
                        time_arr.append(cumulative_time)
                        memory_arr.append(current_max_memory)

                    print(time_arr)
                    print(DTM + " " + dataset + "-" + subset + "-" + method + "-" + dist + " fn/fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])) + "/" + str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
                    print(DTM + " " + dataset + "-" + subset + "-" + method + "-" + dist + " time: ", np.mean(time_arr))
                    print(DTM + " " + dataset + "-" + subset + "-" + method + "-" + dist + " memory: ", np.mean(memory_arr))
                    data_arr.append(current_data_arr)
                    labels.append(dataset + " DTM " + " " + subset + " " + method + " " + dist)
        meta_arr.append((data_arr, labels))
    return meta_arr

def get_RNASim1000_1000M1_HF_starting():
    dataset_arr = ["RNASim1000", "1000M1_HF"]
    for dataset in dataset_arr:
        for method in ["fasttree", "iqtree2"]:
            cumulative_fn_rate = 0.0
            cumulative_fp_rate = 0.0
            cumulative_memory = 0.0
            cumulative_time = 0.0
            time_arr = []
            for replicate in REPLICATE_MAP[dataset]:
                current_model_tree = MODEL_TREE_MAP[dataset].replace("REPLICATE", replicate)
                current_tree = None
                if dataset == "RNASim1000":
                    current_tree = STARTING_RNASim1000_MAP[method].replace("REPLICATE", replicate).replace("SUBSET", "120")
                elif dataset == "1000S4":
                    current_tree = STARTING_1000S4_MAP[method].replace("REPLICATE", replicate).replace("SUBSET", "120")
                elif dataset == "1000M1_HF":
                    current_tree = STARTING_1000M1_HF_MAP[method].replace("REPLICATE", replicate).replace("SUBSET", "120")
                if not Path(current_tree).is_file():
                    # missing_count += 1
                    print(current_tree + " is missing")
                    continue
                nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree, current_tree)
                cumulative_fn_rate += (fn / ei1)
                cumulative_fp_rate += (fp / ei2)

                current_time = None
                current_memory = None
                if(method == "fasttree"):
                    if(dataset == "1000M1_HF"):
                        print(current_tree.split("output")[0] + "errors/fasttree.err")
                        current_time = get_time_fasttree(current_tree.split("output")[0] + "errors/fasttree.err")
                        current_memory = 0
                    else:
                        current_time,current_memory = get_time_memory_from_file(current_tree.split("output")[0] + "errors/fasttree.err")
                elif(method == "iqtree2"):
                    if(dataset == "1000M1_HF"):
                        # print(current_tree.split(".treefile")[0] + ".log")
                        current_time = get_time_iqtree(current_tree.split(".treefile")[0] + ".log")
                        current_memory = get_memory_iqtree(current_tree.split(".treefile")[0] + ".log")
                    else:
                        current_time = get_time_iqtree(current_tree.split("output")[0] + "output/iqtree-full.log")
                        _,current_memory = get_time_memory_from_file(current_tree.split("output")[0] + "errors/iqtree2-complete.err")
                if(method == "mafft"):
                    current_time,current_memory = get_time_memory_from_file(current_tree.split("output")[0] + "errors/mafft.err")
                cumulative_time += current_time
                cumulative_memory += current_memory
                time_arr.append(current_time)

            print("starting " + dataset + "-" + method + " fn/fp rate: ", str(cumulative_fn_rate / len(REPLICATE_MAP[dataset])) + "/" + str(cumulative_fp_rate / len(REPLICATE_MAP[dataset])))
            print("starting " + dataset + "-" + method + " time: ", np.mean(time_arr))
            print("starting " + dataset + "-" + method + " memory: ", str(cumulative_memory /  len(REPLICATE_MAP[dataset])))



# get_incomplete_tree_error()
#print("EXPERIMENT 1")
# get_fasttree_and_internal()
# get_gtm_experiment_1()
#print("EXPERIMENT 2")
# get_gtm_experiment_2()
# print("EXPERIMENT 4")
#, get_gtm_experiment_4()
# print("EXPERIMENT 5")
# get_gtm_experiment_5()
# print("GET 1000M1 Base")
# get_1000M1_base()
#get_200M1_results()
# print("1000M1 PARTIAL TREEMERGE RESULTS")
# get_1000M1_treemerge_results()

# print("EXPERIMENT 6 IQTREE STARTING")
# get_gtm_experiment_6_iqtree_starting()

# print("EXPERIMENT 6")
# get_gtm_experiment6()

# print("EXPERIMENT 6 1000M1 85 support")
# get_gtm_experiment_6_1000M1_support_decompose(85)

# print("EXPERIMENT 6 1000M1 95 support")
# get_gtm_experiment_6_1000M1_support_decompose(95)

# print("EXPERIMENT 6 aa5K centroid")
# get_gtm_experiment_6_aa5K_centroid()

# print("EXPERIMENT 6 aa5K_new 85 support")
# get_gtm_experiment_6_aa5K_support_decompose(85)

# print("EXPERIMENT 6 aa5K_new 95 fasttree")
# get_gtm_experiment_6_fasttree(95)

# print("EXPERIMENT 6 aa5K_new_1000 85 support")
# get_gtm_experiment_6_aa5K_1000_support_decompose(85)

# print("EXPERIMENT 6 aa5K_new_1000 FT_GUIDE_ 85 support")
# get_gtm_experiment_6_aa5K_1000_FT_guide_support_decompose(85)

# print("aa5K_new base")
# get_aa5K_new_base()
# print("EXPERIMENT 6 100k_seeded FT 1000 FT")
# get_gtm_experiment_6_100k_FT_1000_FT()

# print("100k_seeded base")
# get_100k_seeded_base()

# print("100k_seeded_FT_heuristic_1000_FT_GTM")
# get_gtm_experiment_6_100k_FT_heuristic_1000_FT()

#print("Paul RandHetCentro10 base")
#get_randhetcentro10_base()

# print("Paul Constrained-INC")
# get_paul_constrained_inc()

# print("Paul GTM")
# get_paul_gtm()

# print("Paul NJMerge2")
# get_paul_njmerge2()

# print("RandHetCentro10 and 1000M1_HF starting trees")
# get_randhetcentro10_1000M1_starting()
# print("RNASim1000 and 1000M1_HF GTM")
# meta_arr = get_DTM_GTM()
# print("RNASIm1000 and 1000M1_HF Constrained-INC")
# meta_arr = get_DTM_constrained_inc()
# print("RandHetCentro10 and 1000M1_HF TreeMerge")
# meta_arr = get_DTM_treemerge()
# DTM_treemerge_create_box_plots(meta_arr)
# DTM_GTM_create_box_plots(meta_arr)
# DTM_constrained_inc_create_box_plots(meta_arr)

# print("Pair Merge 500 Results")
# meta_arr = get_pairmerge_500_results()
# pair_merge_create_box_plots(meta_arr)

# print("1000S4 RNASim Starting")
# get_RNASim1000_1000S4_starting()
# print("1000S4 Constrained-INC")
# get_DTM_time_memory("Constrained-INC")
# print("RNASim1000 and 1000S4 GTM")
# get_DTM_time_memory("GTM")
# print("Constrained-INC")
# get_DTM_time_memory("Constrained-INC")
# print("RNASim1000 GTM")
# get_DTM_time_memory("GTM")
# print("new treemerge")
# get_new_vs_old_treemerge("FT_NJMerge_TreeMerge")
# print("old treemerge")
# get_new_vs_old_treemerge("OldTreeMerge")
# print("new treemerge DTM RNASim1000 and 1000S4")
# get_DTM_time_memory("FT_NJMerge_TreeMerge")
# print("new treemerge FT_NJMerge")
# get_new_treemerge_time_memory("FT_NJMerge")
# print("new treemerge FT_GTM")
# get_new_treemerge_time_memory("FT_GTM")
# print("new treemerge IQ_GTM")
# get_new_treemerge_time_memory("IQ_GTM")
# print("new treemerge FT_NJMerge")
# get_new_treemerge_time_memory("FT_NJMerge")
print("new treemerge IQ_NJMerge")
get_new_treemerge_time_memory("IQ_NJMerge")
# print("new treemerge FT_Constrained-INC")
# get_new_treemerge_time_memory("FT_Constrained-INC")
# print("new treemerge IQ_Constrained-INC")
get_new_treemerge_time_memory("IQ_Constrained-INC")
print("new treemerge IQ_Constrained-INC")


# print("1000M1_HF RNASim1000 Starting")
# get_RNASim1000_1000M1_HF_starting()
