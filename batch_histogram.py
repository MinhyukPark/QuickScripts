import sys

import matplotlib.pyplot as plt
import numpy as np

import draw_support_histogram

PROJECT = "/projects/tallis/minhyuk2/"
CLUSTER_RUNS = PROJECT + "cluster_runs/"
CLUSTER_RUN_MAP = {
    "Paul/TwoCladesHet/": "PaulTwoCladesHet/PipelineBinning/IQTree/MFP/"
}

def batch_histograms():
    for cluster_run in CLUSTER_RUN_MAP:
        dataset_path = CLUSTER_RUNS + CLUSTER_RUN_MAP[cluster_run]
        replicate_count = 11
        for i in range(1, replicate_count):
            current_replicate = "R0" + str(i)
            if(i == 10):
                current_replicate = "R10"
            current_L_support_tree = dataset_path + current_replicate + "/output/sequence_partition_L.treefile"
            current_R_support_tree = dataset_path + current_replicate + "/output/sequence_partition_L.treefile"

            current_output = "./figures/" + (CLUSTER_RUN_MAP[cluster_run] + current_replicate).replace("/","-")
            left_output = current_output + "-L.png"
            right_output = current_output + "-R.png"
            left_supports = draw_support_histogram.draw_histograms_helper(current_L_support_tree, left_output)
            right_supports = draw_support_histogram.draw_histograms_helper(current_L_support_tree, right_output)

if __name__ == "__main__":
    batch_histograms()
