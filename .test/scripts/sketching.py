#!/usr/bin/env python
"""Experimental implementation of different minhashing approaches.
"""
import mmh3
from functools import partial
from sortedcontainers import SortedList

import dinopy as dp

mmhash = partial(mmh3.hash)


def bottom_sketch(seq, k, s, hf=mmhash):
    """Compute the bottom s sketch of the seqeunce using k-mers.

    Arguments:
        seq (bytes): Seqeunce to be sketched.
        k (int): k-mer size
        s (int): Sketch size, i.e. the number of hash values kept
                in the minimums list.
        hf (function): Hash function to be used. Default mmh3.hash
                with a fixed seed.

    Returns:
        SortedList: the s smallest hash values.
    """
    kmers = dp.qgrams(seq, k)
    # initialize the minimal hash value list with the first s values
    first_s_kmers = [kmers.__next__() for _ in range(s)]
    mins = SortedList([hf(kmer) for kmer in first_s_kmers])
    biggest_min = mins[-1]
    # all_hvs = copy.deepcopy(mins)

    # traverse the remaining kmers
    for index, kmer in enumerate(kmers, s+1):
        hv = hf(kmer)

        # all_hvs.add(hv)  # for debugging
        if hv < biggest_min:
            # print("new min at", index)
            # if the new has value is smaller than the biggest minimizer in
            # the list remove the biggest minimizer and add the new hv to
            # the sorted list
            mins.pop()  # remove the biggest minimizer (at the last position)
            mins.add(hv)
            biggest_min = mins[-1]
        # else:
        #     print("  ", index)

    # print("mins:", all_hvs[:7], min(all_hvs)) # for debugging
    # print("maxs:", all_hvs[-7:]) # for debugging
    return mins


def compare_sketches(sketch_a, sketch_b):
    """Estimate jaccard index by comparing number of shared sketch entries.
    """
    sketch_set_a, sketch_set_b = set(sketch_a), set(sketch_b)
    shared = sketch_set_a & sketch_set_b
    combined = sketch_set_a | sketch_set_b
    return len(shared) / len(combined)
