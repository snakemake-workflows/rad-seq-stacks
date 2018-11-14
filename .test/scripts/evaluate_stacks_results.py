# -*- coding: utf-8 -*-
"""Evaluate a stack mapping versus a simulated rage dataset.
"""
import argparse
import yaml
import sys

import dinopy as dp
import pysam
from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment

from collections import Counter, namedtuple
from functools import partial

import sketching as sk


TSVRecord = namedtuple("TSVRecord", ["locus_id", "seq", "nr_parents",
                                     "nr_snps", "snps", "nr_alleles",
                                     "genotypes"])
VCFRecord = namedtuple("VCFRecord", ["seq", "data"])
GTRecord = namedtuple("GTRecord", ["name", "seq_p5", "seq_p7", "mutations",
                                   "id_reads", "dropout"])
GTStats = namedtuple("GTStats", ["nr_muts", "nr_snps", "nr_inserts",
                                 "nr_deletions", "nr_loci_with_snps"])


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
            orientation = "p7"
        else:
            orientation = "p5"
        return ("SNP", (orientation, snp_pos, (base_from, base_to)))
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
    nr_loci_with_snps = 0

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

        if any((mut_type == "SNP" for mut_type, _ in mutations)):
            nr_loci_with_snps += 1

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

    gt_stats = GTStats(nr_muts, nr_snps, nr_inserts, nr_deletions, nr_loci_with_snps)
    return loc_seqs, gt_stats


def join_seqs(seq_p5, seq_p7):
    """Join two sequences like the ones written by stacks."""
    return "".join((seq_p5, "NNNNNN", dp.reverse_complement(seq_p7)))


def get_stacks_data(args):
    """Read in stacks VCF file."""
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


def evaluate_assembly(assembly, gt_data, stacks_data, gt_stats, args):
    """Analyze an assembly of RAGE vs Stacks data.

    Arguments:
        assembly (dict): Mapping of RAGE sequences to RAGE records
            and associated Stacks records.
        gt_data (list of GTRecords): From a RAGE _gt.yaml file.
        stacks_data (list of TSVRecords): From a Stacks _export.tsv file.
        args (argparse.NameSpace): User parameters.
    """

    nr_of_undiscovered_mutations = 0
    nr_of_discovered_mutations = 0
    # write mapping to file
    with open(args.output, "w") as outfile:
        for (gt_name, gt_seq), (gt_locus, stacks_loci) in assembly.items():
            # user output
            print("Handling gt ", gt_name, ":  ", end="")
            # file output
            print(gt_name, file=outfile)
            print("      ", gt_seq, file=outfile)
            for record in stacks_loci:
                print(f"{', '.join([r.chrom for r in record.data]):>6}",
                      record.seq,
                      file=outfile)

            gt_p5_seq = gt_locus.seq_p5
            gt_p7_seq = gt_locus.seq_p7

            if not stacks_loci:
                print("No matching stack locus found")
                nr_of_undiscovered_mutations += 1
            # compute semiglobal alignments of the loci to verify that
            # they actually match
            for stacks_locus in stacks_loci:
                try:
                    stacks_p5_seq, *_, stacks_p7_seq = stacks_locus.seq.split(b"N")
                except ValueError:
                    stacks_p5_seq = stacks_locus.seq
                    stacks_p7_seq = stacks_locus.seq
                all_p5_alns = align.globalms(gt_p5_seq,
                                             stacks_p5_seq.decode(),
                                             1,   # match score
                                             0,   # mismatch panalty
                                             -5,  # gap open penalty
                                             -3,  # gap extend penalty
                                             penalize_end_gaps=(False, False),
                                             one_alignment_only=True,
                                             )
                all_p7_alns = align.globalms(gt_p7_seq,
                                             dp.reverse_complement(
                                                 stacks_p7_seq.decode()),
                                             1,   # match score
                                             0,   # mismatch panalty
                                             -5,  # gap open penalty
                                             -3,  # gap extend penalty
                                             penalize_end_gaps=(False, False),
                                             one_alignment_only=True,
                                             )
                # pick the first reported alignments
                # these are either unique or good enough
                p5_aln = all_p5_alns[0]
                p7_aln = all_p7_alns[0]
                # print(format_alignment(*p5_aln),
                #       format_alignment(*p7_aln))

                if p5_aln[2] + p7_aln[2] >= 140:
                    print(f"Successful match: {p5_aln[2] + p7_aln[2]}")
                    # print(p5_aln)
                    print([(mutation.alleles, mutation.pos) for mutation in stacks_locus.data])
                    print(gt_locus.mutations)

                    # TODO: Evaluate if right SNPs were found
                    #       (mind that Stacks might not call the allele RAGE
                    #       simulated as root as main allele
                    #         -> consider x>y == y>x)
                    nr_of_discovered_mutations += 1
                    # TODO: Evaluate if the right allele frequencies were
                    #       detected by stacks
                else:
                    print(f"MISMATCH with {stacks_locus.data[0].chrom}")
                    # print(format_alignment(*p5_aln),
                    #       format_alignment(*p7_aln))

            print("\n")
            print("\n", file=outfile)

    # Check to how many RAGE loci each Stacks loci was assigned to.
    # This should always be 1.
    # A value >1 would suggest that a stacks locus is similar to two or
    # more RAGE loci, which is highly unlikely
    locus_assignment_counter = Counter()
    for (gt_name, gt_seq), (gt_locus, stacks_loci) in assembly.items():
        for stacks_locus in stacks_loci:
            # all gt_loci assigned to a stacks locus (identified by
            # the chrom field in the VCF record) should point to
            # the same locus, hence they all just need to be counted as one
            locus_assignment_counter[stacks_locus.data[0].chrom] += 1
    most_assigned_gt_loci = locus_assignment_counter.most_common(5)
    stacks_locus_name, assigned_gt_loci = most_assigned_gt_loci[0]
    if assigned_gt_loci > 1:
        print(f"The 5 most assigned stacks locus id. This should always be 1:"
              f"\n  {most_assigned_gt_loci}\n")

    # find stacks loci that are not in the RAGE loci (singletons and HRLs?)
    # create one entry for each locus assembled by stacks.
    # then find out in the assembly which ones were never assigned to a
    # rage locus

    # 1. create a dictionary with zero for all stacks locus names
    # 2. iterate thourgh the gt_locus -> stacks-locus mapping
    # 3. increment the counter for each stacks locus
    stacks_locus_occurence_count = Counter(
        {rec.data[0].chrom: 0 for rec in stacks_data})
    for (gt_name, gt_seq), (gt_locus, stacks_loci) in assembly.items():
        for locus in stacks_loci:
            # count every gt locus that has an assigned stacks locus
            stacks_locus_occurence_count[locus.data[0].chrom] += 1

    stacks_only_loci = [
        (l, c) for l, c in stacks_locus_occurence_count.items() if c == 0]
    if stacks_only_loci:
        print(f"The following {len(stacks_only_loci)} loci were not simulated "
              "by rage, but identified by stacks. These might include "
              "incompletely digested reads, Null Alleles, Singletons, and "
              "HRLs/ Lumberjack stacks.")
        print([name for name, _ in stacks_only_loci])

    print("")
    print(f"{nr_of_discovered_mutations} loci with mutations were successfully discovered")
    print(f"{nr_of_undiscovered_mutations} loci with mutations were not discovered by stacks")
    print("SNP discovery ratio")
    print(f"{nr_of_discovered_mutations}/{gt_stats.nr_loci_with_snps}")


def main(args):
    """Compare groudn truth with stacks assembly."""
    print(f"Loading gt data", file=sys.stderr)
    gt_data, gt_stats = parse_rage_gt_file(args)

    print(f"Loading stacks data", file=sys.stderr)
    stacks_data = get_stacks_data(args)

    print("Analyzing:", file=sys.stderr)
    assembly = find_matching_loci(gt_data, stacks_data,
                                  similarity=args.similarity_threshold,
                                  )

    print("\n\nLocus Analysis:\n", file=sys.stderr)
    evaluate_assembly(assembly, gt_data, stacks_data, gt_stats, args)

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
