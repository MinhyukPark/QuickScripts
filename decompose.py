from collections import Counter
from multiprocessing import current_process,Pool
import random
import heapq
import sys

import click
import dendropy
from Bio import SeqIO
import scipy
from skbio import Sequence
from skbio.sequence.distance import hamming

import decomposer


# argv guide_tree, subset_size, alignment_on_full_taxa, output_file_prefix

POOL_SIZE=4

def get_closest_sequence_mapping_helper(args):
    return get_closest_sequence_mapping(*args)

def get_closest_sequence_mapping(my_chunk, full_length_sequence_file):
    partial_fragmentary_mapping = {}
    heap_queue = []
    stride_size = 300
    sys.stderr.write("doing work on " + str(len(my_chunk)) + "\n")
    full_length_dict = SeqIO.to_dict(SeqIO.parse(full_length_sequence_file, "fasta"))
    for fragmentary_sequence in my_chunk:
        alignment_sequence_length = len(fragmentary_sequence.seq)
        fragmentary_id = fragmentary_sequence.id
        best_distance = len(fragmentary_sequence)
        best_id = None
        for full_length_sequence_id,full_length_sequence in full_length_dict.items():
            if(best_id is None):
                current_score = scipy.spatial.distance.hamming(list(str(fragmentary_sequence.seq[:stride_size])), list(str(full_length_sequence.seq[:stride_size])))
                to_digest = alignment_sequence_length
                heapq.heappush(heap_queue, (current_score, to_digest, full_length_sequence_id))

        while heapq:
            previous_score,to_digest,current_sequence_id = heapq.heappop(heap_queue)
            if(to_digest < 0):
                best_id = current_sequence_id
                break
            else:
                current_start_index = alignment_sequence_length - to_digest
                current_score = scipy.spatial.distance.hamming(list(str(fragmentary_sequence.seq[current_start_index:current_start_index+stride_size])), list(str(full_length_dict[current_sequence_id].seq[current_start_index:current_start_index+stride_size])))
                heapq.heappush(heap_queue, (previous_score + current_score, to_digest - stride_size, current_sequence_id))
        if(best_id not in partial_fragmentary_mapping):
            partial_fragmentary_mapping[best_id] = []
        partial_fragmentary_mapping[best_id].append(fragmentary_id)
    return partial_fragmentary_mapping


@click.command()
@click.option("--input-tree", required=True, type=click.Path(exists=True), help="The input tree file in newick format")
@click.option("--sequence-file", required=True, type=click.Path(exists=True), help="Aligned sequence file on the full taxa")
@click.option("--fragmentary-sequence-file", required=False, type=click.Path(exists=True), help="Fragmentary aligned sequence file on some set of taxa")
@click.option("--output-prefix", required=True, type=str, help="Output file prefix for each subset")
@click.option("--maximum-size", required=True, type=int, help="Maximum size of output subsets")
@click.option("--longest-edge", is_flag=True, help="Specifying longest edge decomposition")
@click.option("--full-length", required=False, type=int, help="Specifyng a length for the sequences")
@click.option("--longest", required=False, is_flag=True, help="Samples from the longest sequences in the case of incomplete sampling")
@click.option("--incomplete", required=False, type=int, help="Specifying the size of incomplete decomposition (incomplete < maximum-size)")
@click.option("--support-threshold", required=False, type=float, help="Specifying a support threshold for decomposition", default=0.95)
@click.option("--mode", type=click.Choice(["centroid", "support", "heuristic"], case_sensitive=False), required=False, default="centroid")
def decompose_tree(input_tree, sequence_file, fragmentary_sequence_file, output_prefix, maximum_size, longest_edge, full_length, longest, incomplete, support_threshold, mode):
    '''This script decomposes the input tree and outputs induced alignments on the subsets.
    '''
    guide_tree = dendropy.Tree.get(path=input_tree, schema="newick")
    namespace = guide_tree.taxon_namespace
    guide_tree.is_rooted = False
    guide_tree.resolve_polytomies(limit=2)
    guide_tree.collapse_basal_bifurcation()
    guide_tree.update_bipartitions()

    # create mapping of closest sequencs
    fragmentary_mapping = {}
    if(fragmentary_sequence_file != None):
        worker_pool = Pool(processes=POOL_SIZE)
        worker_args_arr = []
        total_fragmentary_sequences = list(SeqIO.parse(open(fragmentary_sequence_file, "r"), "fasta"))
        chunk_size = int(len(total_fragmentary_sequences) / POOL_SIZE)
        for worker_id in range(POOL_SIZE):
            current_chunk = None
            if(chunk_size > 0):
                if(worker_id == POOL_SIZE - 1):
                    current_chunk = total_fragmentary_sequences[int(worker_id*chunk_size):]
                else:
                    current_chunk = total_fragmentary_sequences[int(worker_id*chunk_size):int((worker_id+1)*chunk_size)]
            else:
                current_chunk = total_fragmentary_sequences[int(worker_id):int(worker_id+1)]
            sys.stderr.write("worker id: " + str(worker_id) + "is going to workng on length " + str(len(current_chunk)) + "\n")
            worker_args_arr.append((current_chunk, sequence_file))
        partial_fragmentary_mapping_arr = worker_pool.map(get_closest_sequence_mapping_helper, worker_args_arr)
    sys.stderr.write(str(partial_fragmentary_mapping_arr))
    for partial_fragmentary_mapping in partial_fragmentary_mapping_arr:
        for best_id,fragment_id_arr in partial_fragmentary_mapping.items():
            if (best_id not in fragmentary_mapping):
                fragmentary_mapping[best_id] = []
            fragmentary_mapping[best_id].extend(fragment_id_arr)
    sys.stderr.write(str(fragmentary_mapping))

    trees = None
    sys.stderr.write(str(fragmentary_mapping) + "\n")
    if(fragmentary_sequence_file == None):
        trees = decomposer.decomposeTree(guide_tree, maximum_size, support_threshold, mode=mode)
    elif(fragmentary_sequence_file != None):
        trees = decomposer.fragmentary_decompose_tree(guide_tree, maximum_size, fragmentary_mapping)
    clusters = []
    for tree in trees:
        keep = [n.taxon.label.replace("_"," ") for n in tree.leaf_nodes()]
        clusters.append(set(keep))
    print(len(clusters))

    files = [output_prefix + str(i) + ".out" for i in range(len(clusters))]
    sequence_partitions = [[] for _ in range(len(clusters))]

    for sequence in SeqIO.parse(open(sequence_file), "fasta"):
        for cluster_index,cluster in enumerate(clusters):
            if(sequence.id.replace("_"," ") in cluster):
                sequence_partitions[cluster_index].append(sequence)
                if(sequence.id in fragmentary_mapping):
                    fragmentary_sequence_ids = fragmentary_mapping[sequence.id]
                    for fragmentary_sequence_id in fragmentary_sequence_ids:
                        for fragmentary_sequence in SeqIO.parse(open(fragmentary_sequence_file), "fasta"):
                            if(fragmentary_sequence.id == fragmentary_sequence_id):
                                sequence_partitions[cluster_index].append(fragmentary_sequence)


    for sequence_partition_index,sequence_partition in enumerate(sequence_partitions):
        SeqIO.write(sequence_partition, files[sequence_partition_index], "fasta")


    incomplete_sequences = []
    for sequence_partition_index,sequence_partition in enumerate(sequence_partitions):
        if(incomplete != None):
            if(longest):
                current_list = sorted(sequence_partitions[sequence_partition_index], key=lambda x:len(x.seq.ungap("-")), reverse=True)
                incomplete_sequences.extend(current_list[:incomplete])
            else:
                random.shuffle(sequence_partitions[sequence_partition_index])
                incomplete_sequences.extend(sequence_partitions[sequence_partition_index][:incomplete])
    if(incomplete != None):
        SeqIO.write(incomplete_sequences, output_prefix + "incomplete.out", "fasta")


if __name__ == "__main__":
    decompose_tree()
