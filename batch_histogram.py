import sys

import matplotlib.pyplot as plt
import numpy as np

import draw_support_histogram
import draw_branch_histogram
import scatter_support_branch_length

PROJECT = "/projects/tallis/minhyuk2/"
CLUSTER_RUNS = PROJECT + "cluster_runs/"
CLUSTER_RUN_MAP = {
    "Paul/RandHetCentro10/": "PaulRandHetCentro10/PipelineBinning/IQTree2/MFP/",
    "Paul/TwoCladesHet/": "PaulTwoCladesHet/PipelineBinning/IQTree/MFP/",
}

TREE_ARR = ["ArgS_COG0018-upp_alignment_nuc.fasttree.tre", "RplD_COG0088-upp_alignment_nuc.fasttree.tre", "RplP_COG0197-upp_alignment_nuc.fasttree.tre", "RpsE_COG0098-upp_alignment_nuc.fasttree.tre", "Gtp1_COG0012-upp_alignment_nuc.fasttree.tre", "RplE_COG0094-upp_alignment_nuc.fasttree.tre", "RplR_COG0256-upp_alignment_nuc.fasttree.tre", "RpsG_COG0049-upp_alignment_nuc.fasttree.tre", "HisS_COG0124-upp_alignment_nuc.fasttree.tre", "RplF_COG0097-upp_alignment_nuc.fasttree.tre", "RplV_COG0091-upp_alignment_nuc.fasttree.tre", "RpsH_COG0096-upp_alignment_nuc.fasttree.tre", "PheS_COG0016-upp_alignment_nuc.fasttree.tre", "RplK_COG0080-upp_alignment_nuc.fasttree.tre", "RpoA_COG0202-upp_alignment_nuc.fasttree.tre", "RpsI_COG0103-upp_alignment_nuc.fasttree.tre", "RplA_COG0081-upp_alignment_nuc.fasttree.tre", "RplM_COG0102-upp_alignment_nuc.fasttree.tre", "RpsB_COG0052-upp_alignment_nuc.fasttree.tre", "RplB_COG0090-upp_alignment_nuc.fasttree.tre", "RplN_COG0093-upp_alignment_nuc.fasttree.tre", "RpsC_COG0092-upp_alignment_nuc.fasttree.tre", "RplC_COG0087-upp_alignment_nuc.fasttree.tre", "RplO_COG0200-upp_alignment_nuc.fasttree.tre", "RpsD_COG0522-upp_alignment_nuc.fasttree.tre"]

def batch_histograms():
    for tree in TREE_ARR:
        current_tree_path = "/projects/tallis/minhyuk2/input/Paul/fasttrees_for_min/nidhis/fasttrees-fulltaxa/" + tree
        current_output = "./figures/" + tree
        draw_support_histogram.draw_support_histogram_helper(current_tree_path, current_output + "-support.png")
        draw_branch_histogram.draw_branch_histogram_helper(current_tree_path, current_output + "-branch_length.png")
        scatter_support_branch_length.scatter_support_branch_length_helper(current_tree_path, current_output + "-scatter.png")

if __name__ == "__main__":
    batch_histograms()
