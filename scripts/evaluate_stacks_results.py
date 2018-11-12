# -*- coding: utf-8 -*-
"""Evaluate a stack mapping versus a simulated rage dataset.
"""
import argparse
import yaml
import sys

import dinopy as dp
import pysam
from Bio.pairwise2 import align
# from Bio.pairwise2 import format_alignment

from collections import Counter, namedtuple
from functools import partial

import sketching as sk


TSVRecord = namedtuple("TSVRecord", ["locus_id", "seq", "nr_parents",
                                     "nr_snps", "snps", "nr_alleles",
                                     "genotypes"])
VCFRecord = namedtuple("VCFRecord", ["seq", "data"])
GTRecord = namedtuple("GTRecord", ["name", "seq_p5", "seq_p7", "mutations",
                                   "id_reads", "dropout"])


def normalize_mutation(mut):
    """Normalize SNPs and call mutation type.

    Consensus and mutation might not be clear so if the
    simulation was C>T, a call of T>C is considered valid.

    Returns:
        tuple: Type ("SNP" or "Indel") and (base_from, base_to) for SNPs,
        None for indels
    """
    pos_str, mut_str = mut.split(":")
    read, pos = pos_str.split("@")
    if ">" in mut_str:
        base_from, base_to = mut_str.split(">")
        snp_pos = int(pos)
        if read == "p7":
            base_from = dp.complement(base_from)
            base_to = dp.complement(base_to)
            # 81 is the read length, offset to compensate for rev comp p7 reads
            snp_pos = 81 - snp_pos
        return ("SNP", (snp_pos, (base_from, base_to)))
    else:
        return "Indel", None


def parse_rage_gt_file(args):
    """Read in a RAGE ground truth file.

    Returns:
        list: of GTRecord named tuples each of which hash the entries
        'name', 'seq_p5', 'seq_p7', 'mutations', id_reads', 'dropout'
    """
    with open(args.yaml, 'r') as stream:
        try:
            # read all documents in the data
            inds, loci, *_, hrls = list(yaml.load_all(stream))
        except yaml.YAMLError as exc:
            print(exc)
    nr_muts, nr_snps, nr_inserts, nr_deletions = (0, 0, 0, 0)

    loc_seqs = []
    # filter out all loci with only one allele, i.e. all unmutated loci
    loci_with_snps = ((n, l) for (n, l) in loci.items()
                      if len(l["allele coverages"]) > 1)
    for name, locus in loci_with_snps:
        dropout = []
        mutations = set()
        for n, ind in locus["individuals"].items():
            if ind:
                # dropout events get an empty dict,
                # hence everything that does not evaluate to False
                # is a valid entry with one or two alleles
                for _, allele in ind.items():
                    if allele["mutations"]:
                        normalized_mutations = set(normalize_mutation(a) for a in allele["mutations"])
                        mutations |= normalized_mutations  # extend set
                dropout.append(False)
            else:
                dropout.append(True)

        # compile and append a record for this locus
        id_reads = locus["id reads"]
        # gt_record = GTRecord(name, seq, mutations, id_reads, dropout)
        gt_record = GTRecord(name, locus["p5 seq"], locus["p7 seq"],
                             mutations, id_reads, dropout)
        nr_muts += len(mutations)
        nr_snps += len([mut for mut in mutations if ">" in mut])
        nr_inserts += len([mut for mut in mutations if "+" in mut])
        nr_deletions += len([mut for mut in mutations if "-" in mut])
        loc_seqs.append(gt_record)

    if args.verbose:
        print("Nr of added muts:", nr_muts)
        print("Nr of added snps:", nr_snps)
        print("Nr of added inserts:", nr_inserts)
        print("Nr of added deletions:", nr_deletions)

    return loc_seqs


def join_seqs(seq_p5, seq_p7):
    """Join two sequences like the ones written by stacks."""
    return "".join((seq_p5, "NNNNNN", dp.reverse_complement(seq_p7)))


def get_stacks_data(args):
    """Read in stacks VCF file.
    """
    loc_seqs = []
    haplotypes_file = pysam.VariantFile(args.stacks_haplo, 'r')
    indexed_far = dp.FastaReader(args.stacks_fa)

    record = None
    last_locus = None
    # merge consecutive lines describing SNPs on the same locus.
    for variant_record in haplotypes_file:
        chromosome = variant_record.chrom
        seq = list(indexed_far[chromosome])[0].sequence

        if record is None:
            record = VCFRecord(seq, [variant_record])
            last_locus = variant_record.chrom
        elif variant_record.chrom == last_locus:
            record.data.append(variant_record)
        else:
            loc_seqs.append(record)
            record = VCFRecord(seq, [variant_record])
            last_locus = variant_record.chrom

    return loc_seqs


def find_matching_loci(gt_data, stacks_data, similarity, verbose=True):
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
    k = 6
    s = 50
    sketch = partial(sk.bottom_sketch, k=k, s=s)

    # initialize assembly
    assembly = {
        (gt_record.name, join_seqs(gt_record.seq_p5, gt_record.seq_p7)):
        (gt_record, []) for gt_record in gt_data}
    # print("Sketching stacks data")
    # pre-sketch the data to avoid quadratic (re-) sketching
    sketched_stacks_data = [(sketch(stacks_record.seq), stacks_record)
                            for stacks_record in stacks_data]
    # print("Searching")
    for index, (gt_record) in enumerate(gt_data):
        # print user output
        if verbose and index % 100 == 0:
            print(index)
        joined_gt_seq = join_seqs(gt_record.seq_p5, gt_record.seq_p7)
        s_gt = sketch(joined_gt_seq)
        # compare all stacks locus sketches against the active RAGE loc. sketch
        # Append the loci of all similar sequences to the assembly list.
        for sketch_stacks, record in sketched_stacks_data:
            if sk.compare_sketches(s_gt, sketch_stacks) > similarity:
                assembly[(gt_record.name, joined_gt_seq)][1].append(record)

    return assembly


def evaluate_assembly(assembly, gt_data, stacks_data, args):
    """Analyze an assembly of RAGE vs Stacks data.

    Arguments:
        assembly (dict): Mapping of RAGE sequences to RAGE records
            and associated Stacks records.
        gt_data (list of GTRecords): From a RAGE _gt.yaml file.
        stacks_data (list of TSVRecords): From a Stacks _export.tsv file.
        args (argparse.NameSpace): User parameters.
    """

    # write mapping to file
    if args.output:
        with open(args.output, "w") as outfile:
            for (gt_name, gt_seq), (gt_locus, stacks_loci) in assembly.items():
                print("handling", gt_name)
                print(gt_name, file=outfile)
                print("      ", gt_seq, file=outfile)
                for record in stacks_loci:
                    print("{:>6}".format(", ".join([r.chrom for r in record.data])), record.seq,
                          file=outfile)

                gt_p5_seq = gt_locus.seq_p5
                gt_p7_seq = gt_locus.seq_p7
                # compute semiglobal alignments of the loci to verify that
                # they actually match
                for stacks_locus in stacks_loci:
                    stacks_p5_seq, *_, stacks_p7_seq = stacks_locus.seq.split(b"N")
                    p5_alignments = align.globalms(gt_p5_seq,
                                                   stacks_p5_seq.decode(),
                                                   1,   # match score
                                                   0,   # mismatch panalty
                                                   -5,  # gap open penalty
                                                   -3,  # gap extend penalty
                                                   penalize_end_gaps=(False, False),
                                                   )
                    p7_alignments = align.globalms(gt_p7_seq,
                                                   dp.reverse_complement(stacks_p7_seq.decode()),
                                                   1,   # match score
                                                   0,   # mismatch panalty
                                                   -5,  # gap open penalty
                                                   -3,  # gap extend penalty
                                                   penalize_end_gaps=(False, False),
                                                   )

                    p5_aln = p5_alignments[0]
                    p7_aln = p7_alignments[0]
                    # print(format_alignment(*p5_aln),
                    #       format_alignment(*p7_aln))
                    if p5_aln[4] + p7_aln[4] >= 140:
                        print("Successful match")
                        print([locus.alleles for locus in stacks_locus.data])
                        print(gt_locus.mutations)

                        # TODO: Evaluate if right SNPs were found (mind that Stacks might not call the allele RAGE simulated as root as main allele -> consider x>y == y>x)
                        # TODO: Evaluate if the right allele frequencies were detected by stacks
                print("\n", file=outfile)

    # TODO: check how many mutations were not detected

    # identify length profile of the mapping
    # i.e. how many Stacks loci were assigned to each RAGE locus
    locus_length_counter, locus_id_counter = Counter(), Counter()
    split_loci, id_loci = 0, 0
    total_len = 0
    for (gt_name, gt_seq), (gt_locus, stacks_loci) in assembly.items():
        loc_len = len(stacks_loci)
        total_len += loc_len
        locus_length_counter[loc_len] += 1
        if loc_len > 1:
            split_loci += 1
            if gt_locus.id_reads > 0:
                id_loci += 1
                locus_id_counter[loc_len] += 1

    print(f"Total length: {total_len}")
    print(f"Counts of how many stacks loci were assigned to a RAGE locus:"
          f"\n  {sorted(locus_length_counter.items())}")
    print(f"Of {split_loci} split up loci, {id_loci} had ID reads."
          f"\n  {locus_id_counter}\n")

    # Check to how many RAGE loci each Stacks loci was assigned.
    # This should always be 1.
    # A value >1 would suggest that a stacks locus is similar to two or
    # more RAGE loci, which is highly unlikely
    locus_assignment_counter = Counter()
    for (gt_name, gt_seq), (gt_locus, stacks_loci) in assembly.items():
        for locus in stacks_loci:
            locus_assignment_counter[locus.data.chrom] += 1
    print(f"The 5 most assigned stacks locus id. This should always be 1:"
          f"\n  {locus_assignment_counter.most_common(5)}\n")

    # find stacks loci that are not in the RAGE loci (singletons and HRLs?)
    occurence_count = Counter({rec.data.chrom: 0 for rec in stacks_data})
    for (gt_name, gt_seq), (gt_locus, stacks_loci) in assembly.items():
        for locus in stacks_loci:
            occurence_count[locus.data.chrom] += 1

    # Assemble previously unassigned loci
    loci_to_postprocess = []
    for record in stacks_data:
        if occurence_count[record.data.chrom] == 0:
            loci_to_postprocess.append(record)
    nr_unassigned = sum([1 if count == 0 else 0 for _, count
                         in occurence_count.items()])
    print(f"Number of stacks loci that were not associated with a rage locus: "
          f"{nr_unassigned}")

    # create a secondary assembly with a lower similarity threshold.
    secondary_sim = 0.2
    secondary_assembly = find_matching_loci(gt_data,
                                            loci_to_postprocess,
                                            similarity=secondary_sim,
                                            verbose=False)
    kind_of_similar = 0
    for (gt_name, gt_seq), (gt_locus, stacks_loci) in secondary_assembly.items():
        if stacks_loci:
            print(stacks_loci)
            kind_of_similar += 1
    print(f"Of the {nr_unassigned} unassigned loci in the first pass "
          f"(similarity = {args.similarity_threshold}), {kind_of_similar} "
          f"could be assigned with similarity = {secondary_sim}")

    # Write secondary assembly to file
    if args.output:
        with open(args.output, "a") as outfile:
            print("\n"*20, file=outfile)
            print("Secondary assembly with similarity", file=outfile)
            print("\n", file=outfile)
            for (gt_name, gt_seq), (gt_locus, stacks_loci) in secondary_assembly.items():
                if stacks_loci:
                    print(gt_name, file=outfile)
                    print("      ", gt_seq, file=outfile)
                    for record in stacks_loci:
                        print("{:>6}".format(record.locus_id), record.seq,
                              file=outfile)
                    print("\n", file=outfile)


def main(args):
    """Compare groudn truth with stacks assembly."""
    print(f"Loading gt data", file=sys.stderr)
    gt_data = parse_rage_gt_file(args)

    print(f"Loading stacks data", file=sys.stderr)
    stacks_data = get_stacks_data(args)

    print("Analyzing:", file=sys.stderr)
    assembly = find_matching_loci(gt_data, stacks_data,
                                  similarity=args.similarity_threshold)

    print("\n\nLocus Analysis:\n", file=sys.stderr)
    evaluate_assembly(assembly, gt_data, stacks_data, args)

    # print("\n\nSNPs Analysis:\n", file=sys.stderr)
    # evaluate_snps(assembly, gt_data, stacks_data, args)


def get_argparser():
    """Manage user parameters"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--output",
        help="If an output file should be written",
        default=False,
        )
    parser.add_argument(
        "-g", "-y", "--ground-truth", "--yaml",
        help="Path the a YAML gt file",
        default="RAGEdataset_ATCACG_gt.yaml",
        dest="yaml",
        )
    parser.add_argument(
        "-s", "--stacks-snps-file",
        help="Path to a stacks snps vcf file",
        dest="stacks_haplo",
        )
    parser.add_argument(
        "-f", "--stacks-fasta-file",
        help="Path to a stacks catalog fasta file",
        dest="stacks_fa",
        )
    parser.add_argument(
        "-v", "--verbose",
        help="Print supplementary information",
        dest="verbose",
        action="store_true",
        default=False,
        )
    parser.add_argument(
        "-t", "--similarity-threshold",
        help="The minimal estimated Jaccard similarity for two loci to be"
             "considered similar.",
        dest="similarity_threshold",
        type=float,
        default=0.4,
    )
    return parser


if __name__ == '__main__':
    parser = get_argparser()
    args = parser.parse_args()
    main(args)
