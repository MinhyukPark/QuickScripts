'''
Created on Oct 7, 2019

@author: Vlad
'''


import numpy as np
import dendropy
from dendropy.utility import bitprocessing

import copy
import os
import sys

def decomposeTree(tree, max_subset_size, support_threshold, mode):
    num_leaves = len(tree.leaf_nodes())
    return decompose_tree_helper(num_leaves, tree, max_subset_size, support_threshold, mode)

def decompose_tree_helper(num_taxa, tree, max_subset_size, support_threshold, mode):
    numLeaves = len(tree.leaf_nodes())
    if mode == "heuristic" and numLeaves > max_subset_size * 1.5:
        e = getBestHeuristicEdge(tree, max_subset_size, num_taxa)
        t1, t2 = bipartitionByEdge(tree, e)
        return decompose_tree_helper(num_taxa, t1, max_subset_size, support_threshold, mode) + decompose_tree_helper(num_taxa, t2, max_subset_size, support_threshold, mode)
    elif mode != "heuristic" and numLeaves > max_subset_size:
        if mode == "centroid":
            e = getCentroidEdge(tree)
        elif mode == "random":
            e = getCentroidEdgeRandom(tree, max_subset_size/3)
        elif mode == "longest":
            e = getLongestEdge(tree)
        elif mode == "support":
            e = getSupportEdge(tree, support_threshold)

        t1, t2 = bipartitionByEdge(tree, e)
        return decompose_tree_helper(num_taxa, t1, max_subset_size, support_threshold, mode) + decompose_tree_helper(num_taxa, t2, max_subset_size, support_threshold, mode)
    else:
        if numLeaves >= 1:
            return [tree]
        else:
            sys.exit("tree has fewer than 1 leaves!")

def fragmentary_decompose_tree(tree, max_subset_size, fragmentary_mapping):
    num_full_length_leaves = len(tree.leaf_nodes())
    return fragmentary_decompose_tree_helper(num_full_length_leaves, tree, max_subset_size, fragmentary_mapping)

def fragmentary_decompose_tree_helper(num_taxa, tree, max_subset_size, fragmentary_mapping):
    '''This recursive helper function will return an array of trees
    if the subset size stopping criterion is reached it returns a [tree]
    otherwise it will call 2 more instances of itself and append the subresults together
    '''
    num_full_length_leaves = len(tree.leaf_nodes())
    num_fragmentary_leaves = get_num_fragmentary_leaves(tree.leaf_nodes(), fragmentary_mapping)
    num_leaves = num_full_length_leaves + num_fragmentary_leaves
    # sys.stderr.write("num_full_length_leaves: " + str(num_full_length_leaves))
    # sys.stderr.write("num_fragmentary_leaves: " + str(num_fragmentary_leaves))
    # sys.stderr.write("num_leaves: " + str(num_leaves))
    if(num_leaves > max_subset_size):
        e = get_fragmentary_centroid_edge(tree, max_subset_size, fragmentary_mapping)
        t1, t2 = bipartitionByEdge(tree, e)
        return fragmentary_decompose_tree_helper(num_taxa, t1, max_subset_size, fragmentary_mapping) + fragmentary_decompose_tree_helper(num_taxa, t2, max_subset_size, fragmentary_mapping)
    else:
        if num_leaves >= 1:
            sys.stderr.write("num_full_length_leaves: " + str(num_full_length_leaves) + "\n")
            sys.stderr.write("num_fragmentary_leaves: " + str(num_fragmentary_leaves) + "\n")
            sys.stderr.write("num_leaves: " + str(num_leaves) + "\n")
            return [tree]
        else:
            sys.exit("tree has fewer than 1 leaves!")

def get_num_fragmentary_leaves(leaf_nodes, fragmentary_mapping, is_taxon=False):
    num = 0
    for leaf_node in leaf_nodes:
        if(is_taxon):
            if(leaf_node.label in fragmentary_mapping):
                num += len(fragmentary_mapping[leaf_node.label])
        else:
            if(leaf_node.taxon.label in fragmentary_mapping):
                num += len(fragmentary_mapping[leaf_node.taxon.label])
    return num


def get_fragmentary_centroid_edge(tree, max_subset_size, fragmentary_mapping):
    num_full_length_leaves = len(tree.leaf_nodes())
    num_fragmentary_leaves = get_num_fragmentary_leaves(tree.leaf_nodes(), fragmentary_mapping)
    num_leaves = num_full_length_leaves + num_fragmentary_leaves
    best_balance = float('inf')
    best_edge = None
    sys.stderr.write("searching for best edge in num leaves:")
    sys.stderr.write(str(num_leaves) + str("\n"))
    t1_num_leaves = 0
    t2_num_leaves = 0
    edge_counter = 0
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        # candidate_split_tree = copy.deepcopy(tree)
        # t1, t2 = bipartitionByEdge(candidate_split_tree, edge)
        # t1_full_length_leaves = len(t1.leaf_nodes())
        # t1_fragmentary_leaves = get_num_fragmentary_leaves(t1.leaf_nodes(), fragmentary_mapping)
        # t1_num_leaves = t1_full_length_leaves + t1_fragmentary_leaves
        # t2_full_length_leaves = len(t2.leaf_nodes())
        # t2_fragmentary_leaves = get_num_fragmentary_leaves(t2.leaf_nodes(), fragmentary_mapping)
        # t2_num_leaves = t2_full_length_leaves + t2_fragmentary_leaves
        # balance = max(t1_num_leaves, t2_num_leaves) - min(t1_num_leaves, t2_num_leaves)
        t1_full_length_leaves = edge.bipartition.leafset_taxa(tree.taxon_namespace)
        t1_fragmentary_leaves = get_num_fragmentary_leaves(t1_full_length_leaves, fragmentary_mapping, is_taxon=True)
        t1_num_leaves = len(t1_full_length_leaves) + t1_fragmentary_leaves
        t2_full_length_leaves = get_filtered_leaves(t1_full_length_leaves, tree.leaf_nodes())
        t2_fragmentary_leaves = get_num_fragmentary_leaves(t2_full_length_leaves, fragmentary_mapping, is_taxon=False)
        t2_num_leaves = len(t2_full_length_leaves) + t2_fragmentary_leaves
        sys.stderr.write("searchingc edge " + str(edge_counter) + ": ")
        sys.stderr.write(str(len(t1_full_length_leaves)) + "+" + str(t1_fragmentary_leaves) + "=" + str(t1_num_leaves) + "/" + str(len(t2_full_length_leaves)) + "+" + str(t2_fragmentary_leaves) + "=" + str(t2_num_leaves) + str("\n"))
        assert t1_num_leaves + t2_num_leaves == num_leaves, "Num leaves not adding up"
        balance = abs(num_leaves/2 - t1_num_leaves)
        sys.stderr.write("above edge has balance " + str(balance) + "\n")
        t1_balance = abs(num_leaves/2 - t1_num_leaves)
        t2_balance = abs(num_leaves/2 - t2_num_leaves)
        edge_counter += 1

        if balance < best_balance:
            sys.stderr.write("found better balance " + str(balance) + " and edge " + str(edge_counter) + ": ")
            sys.stderr.write(str(t1_num_leaves) + "/" + str(t2_num_leaves) + str("\n"))
            best_balance = balance
            best_edge = edge
    # sys.stderr.write(str(best_edge.head_node))
    # sys.stderr.write(str(best_edge.length))
    # sys.stderr.write(str(best_edge.head_node.label))
    sys.stderr.write("best balance " + str(best_balance) + ": ")
    sys.stderr.write(str(t1_num_leaves) + "/" + str(t2_num_leaves) + str("\n"))
    return best_edge

def get_filtered_leaves(full_length_leaves, leaf_nodes):
    filtered_arr = []
    full_length_leaf_set = set([leaf.label for leaf in full_length_leaves])
    for leaf_node in leaf_nodes:
        if(leaf_node.taxon.label not in full_length_leaf_set):
            filtered_arr.append(leaf_node)
    return filtered_arr


def bipartitionByEdge(tree, edge):
    newRoot = edge.head_node
    edge.tail_node.remove_child(newRoot)
    newTree = dendropy.Tree(seed_node=newRoot, taxon_namespace = tree.taxon_namespace)
    tree.update_bipartitions()
    newTree.update_bipartitions()
    return tree, newTree

def getCentroidEdge(tree):
    numLeaves = bitprocessing.num_set_bits(tree.seed_node.tree_leafset_bitmask)
    # numLeaves = len(tree.seed_node.leaf_nodes())
    bestBalance = float('inf')
    sys.stderr.write("searching for best edge in num leaves:")
    sys.stderr.write(str(numLeaves) + str("\n"))
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        balance = abs(numLeaves/2 - bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask))
        sys.stderr.write("current_balance:")
        sys.stderr.write(str(balance) + "\n")
        if balance < bestBalance:
            bestBalance = balance
            bestEdge = edge
    sys.stderr.write(str(bestEdge.head_node))
    sys.stderr.write(str(bestEdge.length))
    sys.stderr.write(str(bestEdge.head_node.label))
    return bestEdge

def getSupportEdge(tree, support_threshold):
    numLeaves = bitprocessing.num_set_bits(tree.seed_node.tree_leafset_bitmask)
    # numLeaves = len(tree.seed_node.leaf_nodes())
    bestBalance = float('inf')
    sys.stderr.write("searching for best edge in num leaves:")
    sys.stderr.write(str(numLeaves) + str("\n"))
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        balance = abs(numLeaves/2 - bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask))
        sys.stderr.write("current_balance:")
        sys.stderr.write(str(balance) + "\n")
        if balance < bestBalance and edge.head_node.label is not None and float(edge.head_node.label) > support_threshold:
            bestBalance = balance
            bestEdge = edge
    sys.stderr.write(str(bestEdge.head_node) + "\n")
    sys.stderr.write(str(bestEdge.length) + "\n")
    sys.stderr.write(str(bestEdge.head_node.label) + "\n")
    return bestEdge

def getCentroidEdgeRandom(tree, minBound = 5):
    fullMask = tree.seed_node.tree_leafset_bitmask
    numLeaves = bitprocessing.num_set_bits(fullMask)
    candidates = []
    for edge in tree.postorder_internal_edge_iter():
        if edge.tail_node is None:
            continue

        mask = edge.bipartition.leafset_bitmask
        numMask1 = bitprocessing.num_set_bits(mask)
        numMask2 = numLeaves - numMask1

        if numMask1 >= minBound and numMask2 >= minBound:
            candidates.append(edge)
    return np.random.choice(candidates)

def getLongestEdge(tree):
    numLeaves = bitprocessing.num_set_bits(tree.seed_node.tree_leafset_bitmask)
    # numLeaves = len(tree.seed_node.leaf_nodes())
    longeth_edge = None
    longest_edge_length = 0
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        current_length = edge.length
        if longest_edge_length < current_length:
            longest_edge_length = current_length
            longest_edge = edge
    return longest_edge


def getBestHeuristicEdge(tree, max_subset_size, num_taxa):
    num_leaves = bitprocessing.num_set_bits(tree.seed_node.tree_leafset_bitmask)
    current_best_score = -1
    current_best_edge = None
    subset_L_size = 0
    subset_r_size = 0
    best_L_score = 0
    best_R_score = 0
    best_subset_L_size = 0
    best_subset_R_size = 0
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        if edge.head_node.label is not None:
            subset_L_size = bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask)
            subset_R_size = num_leaves - bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask)
            # sys.stderr.write(str(subset_L_size) + ":" + str(subset_R_size) + "\n")
            current_L_score = heuristic(float(edge.head_node.label), subset_L_size, max_subset_size, num_taxa)
            current_R_score = heuristic(float(edge.head_node.label), subset_R_size, max_subset_size, num_taxa)
            if current_best_score < sum([current_L_score, current_R_score]):
                best_L_score = current_L_score
                best_R_score = current_R_score
                best_subset_L_size = subset_L_size
                best_subset_R_size = subset_R_size
                current_best_score = sum([current_L_score, current_R_score])
                current_best_edge = edge

    # sys.stderr.write(str(current_best_edge.head_node) + "\n")
    # sys.stderr.write(str(current_best_edge.length) + "\n")
    # sys.stderr.write(str(current_best_edge.head_node.label) + "\n")
    sys.stderr.write(str(best_subset_L_size) + ":" + str(best_subset_R_size) + "\n")
    sys.stderr.write(str(best_L_score) + ":" + str(best_R_score) + "\n")
    return current_best_edge


def heuristic(support_value, subset_size, max_subset_size, num_taxa):
    support_component = None
    size_component = None

    support_component = (support_value)
    size_component = 1 - (abs(subset_size - max_subset_size) / num_taxa)
    raw_value = (0.6 * support_component) + (0.4 * size_component)

    if support_value < 0.6:
        raw_value = 0.01 * support_value
    if subset_size < max_subset_size * 0.75:
        raw_value = -1

    return raw_value
    # return max(0, min(0.9999, raw_value))
