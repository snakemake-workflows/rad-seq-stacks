# -*- coding: utf-8 -*-
"""Evaluate a stack mapping versus a simulated rage dataset.
"""
import argparse

import sys
import vcfpy
from collections import defaultdict

# from Bio.pairwise2 import align
# from Bio.pairwise2 import format_alignment

from collections import Counter
from functools import partial

from yaml import dump
from yaml import CDumper as Dumper

import sketching as sk
import file_parser


def join_seqs(seq_p5, seq_p7, join_seq):
    """Join two sequences like the ones written by stacks."""
    return "".join((seq_p5, join_seq, seq_p7))


def find_matching_loci(gt_data, stacks_data, similarity, join_seq,
                       verbose=True):
    """Compare clusterings with minhash sketching.

    Arguments:
        gt_data (list of GTRecords): From a RAGE _gt.yaml file.
        stacks_data (list of TSVRecords): From a Stacks _export.tsv file.
        similarity (float): Similarity threshold for the identification of
            similar sequences. Default: 0.3

    Returns:
        dict: Mapping a tuple of (RAGE locus name, RAGE locus reference
        sequence) to the RAGE locus record and a list of associated Stacks
        locus records.
    """
    # sketching parameters
    k = 7
    s = 100
    sketch = partial(sk.bottom_sketch, k=k, s=s)

    # initialize assembly
    assembly = {
        (gt_record.name, join_seqs(gt_record.seq_p5, gt_record.seq_p7,
                                   join_seq)):
        (gt_record, []) for gt_record in gt_data}
    # print("Sketching stacks data")
    # pre-sketch the data to avoid quadratic (re-) sketching
    sketched_stacks_data = [(sketch(stacks_record.seq), stacks_record)
                            for stacks_record in stacks_data]

    name_mapping = defaultdict(list)

    for index, (gt_record) in enumerate(gt_data):
        # print user output
        if verbose and index % 100 == 0:
            print(index, file=sys.stderr)
        joined_gt_seq = join_seqs(gt_record.seq_p5, gt_record.seq_p7, join_seq)
        # print(joined_gt_seq)
        s_gt = sketch(joined_gt_seq)

        # compare all stacks locus sketches against the active RAGE loc. sketch
        # Append the loci of all similar sequences to the assembly list.

        for sketch_stacks, record in sketched_stacks_data:
            if sk.compare_sketches(s_gt, sketch_stacks) > similarity:
                # assembly[(gt_record.name, joined_gt_seq)][1].append(record)
                name_mapping[record.name].append(gt_record.name)
                record.found = True

    for _, record in sketched_stacks_data:
        if record.found is False:
            print("Not found:", record.seq)

    return name_mapping


def sort_records(rec):
    if rec.CHROM.isdecimal():
        return int(rec.CHROM), rec.POS
    else:
        return rec.CHROM, rec.POS


def rename_vcf_entries(name_mapping, args):

    reader = vcfpy.Reader.from_path(args.stacks_haplo)
    writer = vcfpy.Writer.from_path(args.renamed_vcf, reader.header)
    renamed_records = []
    undetected_counter = 0

    for record in reader:
        stacks_locus = record.CHROM
        try:
            record.CHROM = name_mapping[stacks_locus][0].split(" ")[1]
        except KeyError:
            record.CHROM = f"Unmatched {undetected_counter}"
            undetected_counter += 1
        renamed_records.append(record)

    for record in sorted(renamed_records, key=sort_records):
        writer.write_record(record)

        # TODO find a mathcing stacks record or rename as unmatched #1
        # Replace chrom with the one from the mapping
        # write back to file
        # If the resulting file is not ordered: write into list,
        # sort list by chrom, write back as bulk



def main(args):
    print(f"Loading ddrage gt data", file=sys.stderr)
    gt_data, gt_stats = file_parser.parse_rage_gt_file(args)

    print(f"Loading stacks vcf", file=sys.stderr)
    stacks_data = file_parser.get_stacks_data(args)

    print("Analyzing:", file=sys.stderr)
    name_mapping = find_matching_loci(gt_data, stacks_data,
                                      similarity=args.similarity_threshold,
                                      join_seq=args.join_seq,
    )

    print("\n\nWriting output:\n", file=sys.stderr)
    rename_vcf_entries(name_mapping, args)

    # print("\n\nSNPs Analysis:\n", file=sys.stderr)
    # evaluate_snps(assembly, gt_data, stacks_data, args)


def get_argparser():
    """Manage user parameters"""
    parser = argparse.ArgumentParser()
    # input
    parser.add_argument(
        help="Path to a YAML gt file",
        default="RAGEdataset_ATCACG_gt.yaml",
        dest="yaml",
        )
    parser.add_argument(
        help="Path to a stacks snps vcf file",
        dest="stacks_haplo",
        )
    parser.add_argument(
        help="Path to a stacks catalog fasta file",
        dest="stacks_fa",
        )
    parser.add_argument(
        "-r", "--read-length",
        help="Total simulated read length of one mate.",
        type=int,
        dest="read_length",
        default=100,
        )
    parser.add_argument(
        "-j", "--join-seq",
        help="Sequence used to join mates.",
        required=True,
        type=str,
        dest="join_seq",
        )

    parser.add_argument(
        "-o", "--output",
        help="If an output file should be written",
        dest="renamed_vcf",
        required=True,
        )
    parser.add_argument(
        "-v", "--verbose",
        help="Print supplementary information",
        dest="verbose",
        action="store_true",
        default=False,
        )
    # analysis parameters
    parser.add_argument(
        "-t", "--similarity-threshold",
        help="The minimal estimated Jaccard similarity for two loci to be"
             "considered similar.",
        dest="similarity_threshold",
        type=float,
        default=0.2,
    )
    return parser


if __name__ == '__main__':
    parser = get_argparser()
    args = parser.parse_args()
    main(args)
