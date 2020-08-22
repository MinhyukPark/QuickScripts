'''
Created on Oct 7, 2019

@author: Vlad
'''


import numpy as np
import dendropy
from dendropy.utility import bitprocessing

import os
import sys

def decomposeTree(tree, maxSubsetSize, mode):
    numLeaves = len(tree.leaf_nodes())
    if numLeaves > maxSubsetSize:
        if mode == "centroid":
            e = getCentroidEdge(tree)
        elif mode == "random":
            e = getCentroidEdgeRandom(tree, maxSubsetSize/3)

        t1, t2 = bipartitionByEdge(tree, e)
        return decomposeTree(t1, maxSubsetSize, mode) + decomposeTree(t2, maxSubsetSize, mode)
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

def getCentroidEdge(tree):
    numLeaves = bitprocessing.num_set_bits(tree.seed_node.tree_leafset_bitmask)
    # numLeaves = len(tree.seed_node.leaf_nodes())
    bestBalance = float('inf')
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None:
            continue
        balance = abs(numLeaves/2 - bitprocessing.num_set_bits(edge.bipartition.leafset_bitmask))
        if balance < bestBalance:
            bestBalance = balance
            bestEdge = edge
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
