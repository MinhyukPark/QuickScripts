'''
Created on Oct 7, 2019

@author: Vlad
'''


import numpy as np
import dendropy
from dendropy.utility import bitprocessing

import os
import sys

def decomposeTree(tree, maxSubsetSize, support_threshold, mode):
    numLeaves = len(tree.leaf_nodes())
    if numLeaves > maxSubsetSize:
        if mode == "centroid":
            e = getCentroidEdge(tree, support_threshold)
        elif mode == "random":
            e = getCentroidEdgeRandom(tree, maxSubsetSize/3)
        elif mode == "longest":
            e = getLongestEdge(tree)
        elif mode == "heuristic":
            e = getBestHeuristicEdge(tree)

        t1, t2 = bipartitionByEdge(tree, e)
        return decomposeTree(t1, maxSubsetSize, support_threshold, mode) + decomposeTree(t2, maxSubsetSize, support_threshold, mode)
    else:
        if numLeaves >= 1:
            return [tree]
        else:
            sys.exit("tree has fewer than 1 leaves!")

def bipartitionByEdge(tree, edge):
    newRoot = edge.head_node
    edge.tail_node.remove_child(newRoot)
    newTree = dendropy.Tree(seed_node=newRoot, taxon_namespace = tree.taxon_namespace)
    tree.update_bipartitions()
    newTree.update_bipartitions()
    return tree, newTree

def getCentroidEdge(tree, support_threshold):
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
    sys.stderr.write(str(bestEdge.head_node))
    sys.stderr.write(str(bestEdge.length))
    sys.stderr.write(str(bestEdge.head_node.label))
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


def getBestHeuristicEdge(tree, support_threshold):
    num_leaves = bitprocessing.num_set_bits(tree.seed_node.tree_leafset_bitmask)
    current_best_score = 0
    current_best_edge = None
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        if edge.head_node.label is not None:
            subset_L_size = bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask))
            subset_R_size = num_leaves - bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask))
            current_L_score = heuristic(float(edge.head_node.label), subset_L_size)
            current_R_score = heuristic(float(edge.head_node.label), subset_R_size)
            if current_best_score < max(current_L_score, current_R_score):
                current_best_score = max(current_L_score, current_R_score)
                current_best_edge = edge

    sys.stderr.write(str(current_best_edge.head_node))
    sys.stderr.write(str(current_best_edge.length))
    sys.stderr.write(str(current_best_edge.head_node.label))
    return current_best_edge

def heuristic(support_value, subset_size):
    support_component = None
    size_component = None

    support_component = (support_value)
    size_component = 1 - (abs(subset_size - TARGET) / TARGET)
    raw_value = (0.9 * support_component) + (0.1 * size_component)

    if support_value < 0.6:
        raw_value = 0.01 * support_value
    if subset_size < 250:
        raw_value = 0

    return max(0, min(0.9999, raw_value))
