from compare_trees import main as compare_trees


TALLIS = "/projects/tallis/minhyuk2/"
TALLIS_INPUT = TALLIS + "input/"
CLUSTER_RUNS = TALLIS + "cluster_runs/"

CLUSTER_RUN_MAP = {
    "SATe/1000S1/": "SATe1000S1/NaivePipeline/120/FastTree/",
    "SATe/1000M1/": "SATe1000M1/NaivePipeline/120/FastTree/",
    "SATe/1000L1/": "SATe1000L1/NaivePipeline/120/FastTree/",
    "Vlad/RNASim/1000/": "VladRNASim1000/NaivePipeline/120/FastTree/",
    "Vlad/RNASim/1000_C_500/": "VladRNASim1000_C_500/NaivePipeline/120/FastTree/",
}

DATASET_MAP = {
    "SATe/1000S1/": "rose.mt",
    "SATe/1000M1/": "rose.mt",
    "SATe/1000L1/": "rose.mt",
    "Vlad/RNASim/1000/": "true_tree.tre",
    "Vlad/RNASim/1000_C_500/": "true_tree.tre",
}


for dataset in DATASET_MAP:
    dataset_path = TALLIS_INPUT + dataset
    cumulative_rf_fasttree_rate = 0.0
    cumulative_rf_120_rate = 0.0
    cumulative_rf_500_rate = 0.0
    cumulative_rf_iqtree_rate = 0.0

    cumulative_fn_fasttree_rate = 0.0
    cumulative_fn_120_rate = 0.0
    cumulative_fn_500_rate = 0.0
    cumulative_fn_iqtree_rate = 0.0
    replicate_count = 20
    for i in range(1, replicate_count):
        current_replicate = "R" + str(i)
        current_replicate_path = dataset_path + current_replicate + "/"
        current_model_tree_path = current_replicate_path + DATASET_MAP[dataset]

        current_tree_fasttree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + current_replicate + "/" + "output/fasttree.out"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_fasttree_path)
        if(fp is not fn):
            print("error: fp != fn")
        # print(nl, ei1, ei2, fp, fn, rf)
        # print("fp == fn: ", fp == fn)
        cumulative_rf_fasttree_rate += rf
        cumulative_fn_fasttree_rate += (fn / ei1)

        current_tree_120_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + current_replicate + "/" + "output/treemerge.tree"
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_120_path)
        if(fp is not fn):
            print("error: fp != fn")
        # print(nl, ei1, ei2, fp, fn, rf)
        # print("fp == fn: ", fp == fn)
        cumulative_rf_120_rate += rf
        cumulative_fn_120_rate += (fn / ei1)

        current_tree_500_path = current_tree_120_path.replace("120", "500")
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_500_path)
        if(fp is not fn):
            print("error: fp != fn")
        # print(nl, ei1, ei2, fp, fn, rf)
        # print("fp == fn: ", fp == fn)
        cumulative_rf_500_rate += rf
        cumulative_fn_500_rate += (fn / ei1)

        current_tree_iqtree_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[dataset] + current_replicate + "/" + "output/IQtree-ModelFinder-Full.treefile"
        current_tree_iqtree_path = current_tree_iqtree_path.replace("NaivePipeline/120/FastTree/", "IQTree/ModelFinder/")
        nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree_iqtree_path)
        if(fp is not fn):
            print("error: fp != fn")
        # print(nl, ei1, ei2, fp, fn, rf)
        # print("fp == fn: ", fp == fn)
        cumulative_rf_iqtree_rate += rf
        cumulative_fn_iqtree_rate += (fn / ei1)

    average_rf_rate = cumulative_rf_120_rate / replicate_count
    average_fn_rate = cumulative_fn_120_rate / replicate_count
    print("120 dataset: ", dataset, end=" - ")
    print("av rf rate: ", average_rf_rate)
    print("av fn rate: ", average_fn_rate)

    average_rf_rate = cumulative_rf_500_rate / replicate_count
    average_fn_rate = cumulative_fn_500_rate / replicate_count
    print("500 dataset: ", dataset, end=" - ")
    print("av rf rate: ", average_rf_rate)
    print("av fn rate: ", average_fn_rate)

    average_rf_rate = cumulative_rf_fasttree_rate / replicate_count
    average_fn_rate = cumulative_fn_fasttree_rate / replicate_count
    print("FastTree dataset: ", dataset, end=" - ")
    print("av rf rate: ", average_rf_rate)
    print("av fn rate: ", average_fn_rate)

    average_rf_rate = cumulative_rf_iqtree_rate / replicate_count
    average_fn_rate = cumulative_fn_iqtree_rate / replicate_count
    print("IQTree dataset: ", dataset, end=" - ")
    print("av rf rate: ", average_rf_rate)
    print("av fn rate: ", average_fn_rate)
