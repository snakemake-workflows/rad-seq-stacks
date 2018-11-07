#!/usr/bin/env python
"""Experimental implementation of different minhashing approaches.
"""
import os
import sys
import numpy
import mmh3
import copy
import numpy as np
import random
import argparse
from pprint import pprint
from functools import partial
from sortedcontainers import SortedList

import dinopy as dp

mmhash = partial(mmh3.hash)


def bottom_sketch(seq, k, s, hf=mmhash):
    """Compute the bottom s sketch of the seqeunce using k-mers.

    Arguments:
        seq (bytes): Seqeunce to be sketched.
        k (int): k-mer size
        s (int): Sketch size, i.e. the number of hash values kept in the minimums list.
        hf (function): Hash function to be used. Default mmh3.hash with a fixed seed.

    Returns:
        SortedList: the s smallest hash values 
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
            # if the new has value is smaller than the biggest minimizer in the list
            # remove the biggest minimizer and add the new hv to the sorted list
            old = mins.pop()
            mins.add(hv)
            biggest_min = mins[-1]
        # else:
        #     print("  ", index)

    # print("mins:", all_hvs[:7], min(all_hvs)) # for debugging
    # print("maxs:", all_hvs[-7:]) # for debugging
    return mins


def multiple_hf_random_permutation(hf=mmhash):
    ...


def multiple_hf_xor(seq, k, nr_hfs, seed=42, hf=mmhash):
    """Compute minhash sketch with nr_hfs hash functions.

    Arguments:
        seq (bytes): Seqeunce to be sketched.
        k (int): K-mer size.
        nr_hfs (int): Number of hash functions in the sketch.
        hf (function): Hash function to be used. Default: mmh3.hash with fixed seed
    """
    # use seed to make hf generation predicable for testing
    np.random.seed(seed)
    # compute one 32bit xor value for each hash function
    # a hash function will be computed as g_i(x) = hf(x) ^ xor_values[i]
    xor_values = np.random.randint(0, 4294967296, size=nr_hfs, dtype=np.uint32)
    hv_minimums = np.zeros(nr_hfs, dtype=np.int32)

    # print("XORS", xor_values)
    # print("\n")
    # print(hv_minimums)
    for kmer in dp.qgrams(seq, k):
        hv = hf(kmer)
        hvs = np.bitwise_xor(hv, xor_values)
        hv_minimums = np.minimum(hv_minimums, hvs)
        # print(hv_minimums)

    return hv_minimums


def compare_sketches(sketch_a, sketch_b, verbose=False):
    """Estimate jaccard index by comparing number of shared sketch entries.
    """
    sketch_set_a, sketch_set_b = set(sketch_a), set(sketch_b)
    shared = sketch_set_a & sketch_set_b
    combined = sketch_set_a | sketch_set_b
    if verbose:
        print("Sketch A:\n{}\nSketch B:\n{}".format(sorted(sketch_set_a), sorted(sketch_set_b)))
        print("Shared (A ∩ B): {}\n{}".format(len(shared), sorted(shared)))
        print("Combinde (A ∪ B): {}\n{}".format(len(combined), sorted(combined)))
    return len(shared) / len(combined)


def almost_like(seq, seed, threshold):
    """Create a copy that is almost like the input seqeunce, with an error rate of threshold.

    Arguments:
        seq (bytes): The sequence for whoch an altered copy will be created.
        seed (bool): If set, fix the seed to the global variable SEED (default 42).
        threshold (float): Error probability. Default 0.05 via argparse.

    Returns:
        bytes: A slightly altered seqeunce.
    """
    if seed:
        random.seed(SEED)
        np.random.seed(SEED)
    replacements = b"ARNDCQEGHILKMFPSTWYV"
    rand = np.random.random
    output = [c if rand() > threshold else random.choice(replacements) for c in seq]
    # print(" Input: {}".format(seq))
    # print("Output: {}".format(bytes(output)))
    return bytes(output)


def change_seed(seed):
    SEED = seed
    mmhash = partial(mmh3.hash, seed=SEED)


def main(args):
    # read small test dataset, drop the names 
    far = dp.FastaReader(args.input_file)
    seqs = [record.sequence for record in far.chromosomes()]
    # define compact function call to slightly alter a sequence
    almost = partial(almost_like, seed=args.fix_seed, threshold=args.error_rate)

    print("Bottom s sketch")
    unaltered_bs_sketch = bottom_sketch(seqs[0], k=6, s=5)
    altered_bs_sketch = bottom_sketch(almost(seqs[0]), k=6, s=5)
    print("Unaltered", unaltered_bs_sketch, sep="\n")
    print("Altered", altered_bs_sketch, sep="\n")

    print("XOR")
    unaltered_xor_sketch = multiple_hf_xor(seqs[0], k=6, nr_hfs=5)
    altered_xor_sketch = multiple_hf_xor(almost(seqs[0]), k=6, nr_hfs=5)
    print("Unaltered", unaltered_xor_sketch, sep="\n")
    print("Altered", altered_xor_sketch, sep="\n")



def get_argpaser():
    description = "Test environment for various minhashing schemes."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-i", "--input_file",
        help="Input file. Default: uniprot_35.fasta",
        dest="input_file",
        default="uniprot_35.fasta",
        )
    parser.add_argument(
        "-f", "--fix-seed",
        help="Fix all random seeds to 42.",
        action="store_true",
        default=False,
        )
    parser.add_argument(
        "-e", "--error-rate",
        help="Error rate with which amino acids will be changed. Default: 0.05",
        default=0.05,
        type=float,
        dest="error_rate",
        )
    return parser


if __name__ == '__main__':
    parser = get_argpaser()
    args = parser.parse_args()
    main(args)
