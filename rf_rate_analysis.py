from compare_trees import main as compare_trees


PROJECT = "/projects/sciteam/bbaz/minhyuk2/"
PROJECT_INPUT = PROJECT + "input/"
CLUSTER_RUNS = PROJECT + "cluster_runs/"

CLUSTER_RUN_MAP = {
    "Paul/TwoCladesHet/": "PaulTwoCladesHet/IQTree/"
}

DATASET_MAP = {
    "Paul/TwoCladesHet": "true-tree.tre"
}


for dataset in DATASET_MAP:
    dataset_path = PROJECT_INPUT + dataset
    cumulative_rf_fasttree_rate = 0.0
    cumulative_rf_iqtree_jc_rate = 0.0
    cumulative_rf_iqtree_gtrgamma_rate = 0.0

    cumulative_fn_fasttree_rate = 0.0
    cumulative_fn_iqtree_jc_rate = 0.0
    cumulative_fn_iqtree_gtrgamma_rate = 0.0

    replicate_count = 11
    for i in range(1, replicate_count):
        current_replicate = "R0" + str(i)
        if(i == 10):
            current_replicate = "R10"
        current_replicate_path = dataset_path + current_replicate + "/"
        current_model_tree_path = current_replicate_path + DATASET_MAP[dataset]
        current_tree_fasttree_path = current_replicate_path + "fasttree.tre"

        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_fasttree_path)
        if(fp is not fn):
            print("error: fp != fn")
        cumulative_rf_fasttree_rate += rf
        cumulative_fn_fasttree_rate += (fn / ei1)

        current_tree_iqtree_jc_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "JC69/" + current_replicate + "/output/TwoCladesHetIQTree.treefile"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_iqtree_jc_path)
        if(fp is not fn):
            print("error: fp != fn")
        cumulative_rf_iqtree_jc_rate += rf
        cumulative_fn_iqtree_jc_rate += (fn / ei1)

        current_tree_iqtree_gtrgamma_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + "GTRGamma/" + current_replicate + "/output/TwoCladesHetIQTree.treefile"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_iqtree_gtrgamma_path)
        if(fp is not fn):
            print("error: fp != fn")
        cumulative_rf_iqtree_gtrgamma_rate += rf
        cumulative_fn_iqtree_gtrgamma_rate += (fn / ei1)

    average_rf_rate = cumulative_rf_fasttree_rate / (replicate_count - 1)
    average_fn_rate = cumulative_fn_fasttree_rate / (replicate_count - 1)
    print("FastTree dataset: ", dataset, end=" - ")
    print("av rf rate: ", average_rf_rate)
    print("av fn rate: ", average_fn_rate)

    average_rf_rate = cumulative_rf_iqtree_jc_rate / (replicate_count - 1)
    average_fn_rate = cumulative_fn_iqtree_jc_rate / (replicate_count - 1)
    print("IQTree JC69 dataset: ", dataset, end=" - ")
    print("av rf rate: ", average_rf_rate)
    print("av fn rate: ", average_fn_rate)

    average_rf_rate = cumulative_rf_iqtree_gtrgamma_rate / (replicate_count - 1)
    average_fn_rate = cumulative_fn_iqtree_gtrgamma_rate / (replicate_count - 1)
    print("IQTree GTRGamma dataset: ", dataset, end=" - ")
    print("av rf rate: ", average_rf_rate)
    print("av fn rate: ", average_fn_rate)
