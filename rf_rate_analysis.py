import pathlib
from compare_trees import main as compare_trees


# PROJECT = "/projects/sciteam/bbaz/minhyuk2/"
PROJECT = "/projects/tallis/minhyuk2/"
PROJECT_INPUT = PROJECT + "input/"
CLUSTER_RUNS = PROJECT + "cluster_runs/"

CLUSTER_RUN_MAP = {
    "Paul/TwoCladesHet/": "PaulTwoCladesHet/"
}

DATASET_MAP = {
    "Paul/TwoCladesHet/": "true-tree.tre"
}


for dataset in DATASET_MAP:
    dataset_path = PROJECT_INPUT + dataset

    cumulative_fn_fasttree_rate = 0.0
    cumulative_fn_iqtree_jc_rate = 0.0
    cumulative_fn_iqtree_gtrgamma_rate = 0.0
    cumulative_fn_naive_pipeline_2_rate = 0.0
    cumulative_fn_naive_pipeline_120_rate = 0.0
    cumulative_fn_pipeline_binning_50_rate = 0.0
    cumulative_fn_pipeline_binning_60_rate = 0.0
    cumulative_fn_pipeline_binning_70_rate = 0.0
    cumulative_fn_pipeline_binning_80_rate = 0.0
    cumulative_fn_pipeline_binning_90_rate = 0.0

    replicate_count = 11
    for i in range(1, replicate_count):
        current_replicate = "R0" + str(i)
        if(i == 10):
            current_replicate = "R10"
        current_replicate_path = dataset_path + current_replicate + "/"
        current_model_tree_path = current_replicate_path + DATASET_MAP[dataset]
        current_tree_fasttree_path = current_replicate_path + "fasttree.tre"

        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_fasttree_path)
        if(fp != fn):
            print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
            print("error: fp != fn")
        cumulative_fn_fasttree_rate += (fn / ei1)

        current_tree_iqtree_jc_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "IQTree/JC69/" + current_replicate + "/output/TwoCladeHetIQTree.treefile"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_iqtree_jc_path)
        if(fp != fn):
            print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
            print("error: fp != fn")
        cumulative_fn_iqtree_jc_rate += (fn / ei1)

        current_tree_iqtree_gtrgamma_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "IQTree/GTRGamma/" + current_replicate + "/output/TwoCladeHetIQTree.treefile"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_iqtree_gtrgamma_path)
        if(fp != fn):
            print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
            print("error: fp != fn")
        cumulative_fn_iqtree_gtrgamma_rate += (fn / ei1)

        current_naive_pipeline_2_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "NaivePipeline/2Subsets/FastTree/" + current_replicate + "/output/gtm.tree"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_naive_pipeline_2_tree_path)
        if(fp != fn):
           print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
           print("error: fp != fn")
        cumulative_fn_naive_pipeline_2_rate += (fn / ei1)

        current_naive_pipeline_120_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "NaivePipeline/120/FastTree/" + current_replicate + "/output/gtm.tree"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_naive_pipeline_120_tree_path)
        if(fp != fn):
           print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
           print("error: fp != fn")
        cumulative_fn_naive_pipeline_120_rate += (fn / ei1)

        current_pipeline_binning_50_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/gtm-50.tree"
        current_pipeline_binning_50_tree_file = pathlib.Path(current_pipeline_binning_50_tree_path)
        if not current_pipeline_binning_50_tree_file.exists():
            L_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-50-L.out", "r")
            R_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-50-R.out", "r")
            L_compat_status = L_compat_file.read()
            R_compat_status = R_compat_file.read()
            print(L_compat_status)
            print(R_compat_status)
            if(L_compat_status == "0 0 " and R_compat_status == "0 0 "):
                print("binning strategy with threshold 50 results in complete merge at: " + current_replicate )
            current_pipeline_binning_50_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/sequence_partition_50-merged.tree.treefile"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_pipeline_binning_50_tree_path)
        if(fp != fn):
           print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
           print("error: fp != fn")
        cumulative_fn_pipeline_binning_50_rate += (fn / ei1)


        current_pipeline_binning_60_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/gtm-60.tree"
        current_pipeline_binning_60_tree_file = pathlib.Path(current_pipeline_binning_60_tree_path)
        if not current_pipeline_binning_60_tree_file.exists():
            L_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-60-L.out", "r")
            R_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-60-R.out", "r")
            L_compat_status = L_compat_file.read()
            R_compat_status = R_compat_file.read()
            if(L_compat_status == "0 0 " and R_compat_status == "0 0 "):
                print("binning strategy with threshold 60 results in complete merge at: " + current_replicate )
            current_pipeline_binning_60_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/sequence_partition_60-merged.tree.treefile"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_pipeline_binning_60_tree_path)
        if(fp != fn):
           print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
           print("error: fp != fn")
        cumulative_fn_pipeline_binning_60_rate += (fn / ei1)

        current_pipeline_binning_70_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/gtm-70.tree"
        current_pipeline_binning_70_tree_file = pathlib.Path(current_pipeline_binning_70_tree_path)
        if not current_pipeline_binning_70_tree_file.exists():
            L_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-70-L.out", "r")
            R_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-70-R.out", "r")
            L_compat_status = L_compat_file.read()
            R_compat_status = R_compat_file.read()
            if(L_compat_status == "0 0 " and R_compat_status == "0 0 "):
                print("binning strategy with threshold 70 results in complete merge at: " + current_replicate )
            current_pipeline_binning_70_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/sequence_partition_70-merged.tree.treefile"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_pipeline_binning_70_tree_path)
        if(fp != fn):
           print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
           print("error: fp != fn")
        cumulative_fn_pipeline_binning_70_rate += (fn / ei1)

        current_pipeline_binning_80_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/gtm-80.tree"
        current_pipeline_binning_80_tree_file = pathlib.Path(current_pipeline_binning_80_tree_path)
        if not current_pipeline_binning_80_tree_file.exists():
            L_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-80-L.out", "r")
            R_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-80-R.out", "r")
            L_compat_status = L_compat_file.read()
            R_compat_status = R_compat_file.read()
            if(L_compat_status == "0 0 " and R_compat_status == "0 0 "):
                print("binning strategy with threshold 80 results in complete merge at: " + current_replicate )
            current_pipeline_binning_80_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/sequence_partition_80-merged.tree.treefile"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_pipeline_binning_80_tree_path)
        if(fp != fn):
           print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
           print("error: fp != fn")
        cumulative_fn_pipeline_binning_80_rate += (fn / ei1)

        current_pipeline_binning_90_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/gtm-90.tree"
        current_pipeline_binning_90_tree_file = pathlib.Path(current_pipeline_binning_90_tree_path)
        if not current_pipeline_binning_90_tree_file.exists():
            L_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-90-L.out", "r")
            R_compat_file = open(CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/compatsequence_partition_-90-R.out", "r")
            L_compat_status = L_compat_file.read()
            R_compat_status = R_compat_file.read()
            if(L_compat_status == "0 0 " and R_compat_status == "0 0 "):
                print("binning strategy with threshold 90 results in complete merge at: " + current_replicate )
            current_pipeline_binning_90_tree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "PipelineBinning/IQTree/MFP/" + current_replicate + "/output/sequence_partition_90-merged.tree.treefile"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_pipeline_binning_90_tree_path)
        if(fp != fn):
           print(f"fp: {fp}, fn: {fn}, ei1: {ei1}, ei2: {ei2}")
           print("error: fp != fn")
        cumulative_fn_pipeline_binning_90_rate += (fn / ei1)

    average_fn_rate = cumulative_fn_fasttree_rate / (replicate_count - 1)
    print("FastTree dataset: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)

    average_fn_rate = cumulative_fn_iqtree_jc_rate / (replicate_count - 1)
    print("IQTree JC69 dataset: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)

    average_fn_rate = cumulative_fn_iqtree_gtrgamma_rate / (replicate_count - 1)
    print("IQTree GTRGamma dataset: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)

    average_fn_rate = cumulative_fn_naive_pipeline_2_rate / (replicate_count - 1)
    print("Naive Pipeline 2 Subsets: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)

    average_fn_rate = cumulative_fn_naive_pipeline_120_rate / (replicate_count - 1)
    print("Naive Pipeline 120: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)

    average_fn_rate = cumulative_fn_pipeline_binning_50_rate / (replicate_count - 1)
    print("Pipeline binning 50: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)

    average_fn_rate = cumulative_fn_pipeline_binning_60_rate / (replicate_count - 1)
    print("Pipeline binning 60: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)

    average_fn_rate = cumulative_fn_pipeline_binning_70_rate / (replicate_count - 1)
    print("Pipeline binning 70: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)

    average_fn_rate = cumulative_fn_pipeline_binning_80_rate / (replicate_count - 1)
    print("Pipeline binning 80: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)

    average_fn_rate = cumulative_fn_pipeline_binning_90_rate / (replicate_count - 1)
    print("Pipeline binning 90: ", dataset, end=" - ")
    print("av fn rate: ", average_fn_rate)
