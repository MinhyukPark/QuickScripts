import main as compare_trees from ./compare_trees.py


TALLIS = "/projects/tallis/minhyuk2/"
TALLIS_INPUT = TALLIS + "input/"
CLUSTER_RUNS = TALLIS + "cluster_runs/"

CLUSTER_RUN_MAP = {
    "SATe/1000S1/": "SATe1000S1/IQTree/ModelFinder/",
    "SATe/1000M1/": "SATe1000M1/IQTree/ModelFinder/",
    "SATe/1000L1/": "SATe1000L1/IQTree/ModelFinder/",
    "Vlad/RNASim/1000/": "VladRNASim1000/IQTree/ModelFinder/",
    "Vlad/RNASim/1000_C_500/": "VladRNASim1000_C_500/IQTree/ModelFinder/",
}

DATASET_MAP = {
    "SATe/1000S1/": "rose.mt",
    "SATe/1000M1/": "rose.mt",
    "SATe/1000L1/": "rose.mt",
    "Vlad/RNASim/1000/": "true_tree.tre",
    "Vlad/RNASim/1000_C_500/": "true_tree.tre",
}

is_iqtree = True

for dataset,model_tree in dataset_map:
    dataset_path = TALLIS_INPUT + dataset
    cumulative_rf_rate = 0.0
    replicate_count = 20
    for i in range(replicate_count):
        current_replicate = "R" + str(i)
        current_replicate_path = dataset_path + current_replicate + "/"
        current_model_tree_path = current_replicate_path + model_tree
        if is_iqtree:
            current_tree_path = CLUSTER_RUN_MAP[dataset] + "output/IQTree-ModelFinder-Full.treefile"
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree)
            cumulative_rf_rate += rf
        else:
            current_tree_path = CLUSTER_RUN_MAP[dataset] + "output/IQTree-ModelFinder-Full.RAXMLSOMETHING" # TODO: add raxml here
            nl, ei1, ei2, fp, fn, rf = compare_trees(current_model_tree_path, current_tree)
            cumulative_rf_rate += rf
    average_rf_rate = current_rf_rate / replicate_count
    print(dataset, end="")
    print(average_rf_rate)
