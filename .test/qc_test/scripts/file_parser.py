import sys
import dinopy as dp
import pysam
import yaml
from collections import namedtuple
from recordtype import recordtype

TSVRecord = namedtuple("TSVRecord", ["locus_id", "seq", "nr_parents",
                                     "nr_snps", "snps", "nr_alleles",
                                     "genotypes"])
GTRecord = namedtuple("GTRecord", ["name", "seq_p5", "seq_p7", "mutations",
                                   "id_reads", "dropout", "alleles",
                                   "allele_frequencies"])
GTStats = namedtuple("GTStats", ["nr_muts", "nr_snps", "nr_inserts",
                                 "nr_deletions", "nr_loci_with_snps",
                                 "nr_loci_with_muts", "nr_loci"])

VCFRecord = recordtype("VCFRecord", ["name", "seq", "data", "found"])


def normalize_mutation(mut, offset):
    """Normalize SNPs and call mutation type.

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
            # move SNP to compensate for p7 sequence being affixed to p5 seq
            snp_pos = snp_pos + offset
            orientation = "p7"
        else:
            orientation = "p5"
        return ("SNP", (orientation, snp_pos, (base_from, base_to)))
    else:
        return "Indel", None


def process_mutations(mutations):
    """Pack mutations into a dict
    """
    processed_mutations = []

    for mutation in mutations:
        (_, (orientation, pos, (ref, alt))) = mutation
        processed_mutations.append({
            "orientation": orientation,
            "pos": pos,
            "ref": ref,
            "alt": alt,
        })
    return processed_mutations


def parse_rage_gt_file(args):
    """Read in a RAGE ground truth file.

    Returns:
        list: of GTRecord named tuples each of which hash the entries
        'name', 'seq_p5', 'seq_p7', 'mutations', id_reads', 'dropout'
    """
    with open(args.yaml, 'r') as stream:
        try:
            # read all documents in the data
            inds, loci, *other = list(yaml.load_all(stream, Loader=yaml.FullLoader))
        except yaml.YAMLError as exc:
            print(exc)
    nr_muts, nr_snps, nr_inserts, nr_deletions = 0, 0, 0, 0
    nr_loci_with_snps, nr_loci_with_muts = 0, 0

    p5_enz = list(inds["Individual Information"].values())[0]["p5 overhang"]
    loc_seqs = []
    # filter out all loci with only one allele, i.e. all unmutated loci
    loci_with_snps = ((n, l) for (n, l) in loci.items()
                      if len(l["allele coverages"]) > 1)

    # print("inds", inds)
    spacer_lengths = [len(i["p5 spacer"]) for i
                      in inds["Individual Information"].values()]
    # spacer_variance = max(spacer_lengths) - min(spacer_lengths)
    overhang_lengths = [len(i["p5 overhang"]) for i
                        in inds["Individual Information"].values()]
    # overhang_variance = max(overhang_lengths) - min(overhang_lengths)
    offset = args.read_length - (min(spacer_lengths) + min(overhang_lengths)) \
        + len(args.join_seq)

    for name, locus in loci.items():
        dropout = []
        mutations = set()
        gt_alleles = {}
        for n, ind in locus["individuals"].items():
            if ind:
                # dropout events get an empty dict,
                # hence everything that does not evaluate to False
                # is a valid entry with one or two alleles
                for nr_allele, allele in ind.items():
                    normalized_mutations = set()
                    if allele["mutations"]:
                        normalized_mutations = set(
                            normalize_mutation(a, offset)
                            for a in allele["mutations"]
                        )
                        if nr_allele not in gt_alleles:
                            gt_alleles[nr_allele] = {
                                "frequency": locus["allele frequencies"][nr_allele],
                                "mutations": process_mutations(normalized_mutations),
                                "coverage": locus["allele coverages"][nr_allele],
                            }

                        mutations |= normalized_mutations  # extend set
                dropout.append(False)
            else:
                dropout.append(True)

        if any((mut_type == "SNP" for mut_type, _ in mutations)):
            nr_loci_with_snps += 1
            nr_loci_with_muts += 1
        elif mutations:
            nr_loci_with_muts += 1  # locus with indels only

        # compile and append a record for this locus
        id_reads = locus["id reads"]
        # gt_record = GTRecord(name, seq, mutations, id_reads, dropout)
        gt_record = GTRecord(name, p5_enz + locus["p5 seq"], locus["p7 seq"],
                             mutations, id_reads, dropout,
                             gt_alleles, locus["allele frequencies"])
        nr_muts += len(mutations)
        nr_snps += len([mut for mut in mutations if ">" in mut])
        nr_inserts += len([mut for mut in mutations if "+" in mut])
        nr_deletions += len([mut for mut in mutations if "-" in mut])
        loc_seqs.append(gt_record)

    if args.verbose:
        print("Nr of added muts:", nr_muts, file=sys.stderr)
        print("Nr of added snps:", nr_snps, file=sys.stderr)
        print("Nr of added inserts:", nr_inserts, file=sys.stderr)
        print("Nr of added deletions:", nr_deletions, file=sys.stderr)

    gt_stats = GTStats(nr_muts, nr_snps, nr_inserts, nr_deletions,
                       nr_loci_with_snps, nr_loci_with_muts, len(loci))
    return loc_seqs, gt_stats


def get_stacks_data(args):
    """Read in stacks VCF file."""
    loc_seqs = dict()
    haplotypes_file = pysam.VariantFile(args.stacks_haplo, 'r')
    indexed_far = dp.FastaReader(args.stacks_fa)

    record = None
    last_locus = None
    chromosome = None

    # merge consecutive lines describing SNPs on the same locus.
    for variant_record in haplotypes_file:
        chromosome = variant_record.chrom
        if record is None:
            seq = list(indexed_far[chromosome])[0].sequence
            record = VCFRecord(chromosome, seq, [variant_record], False)
            last_locus = variant_record.chrom
        elif variant_record.chrom == last_locus:
            record.data.append(variant_record)
        else:
            loc_seqs[last_locus] = record
            seq = list(indexed_far[chromosome])[0].sequence
            record = VCFRecord(chromosome, seq, [variant_record], False)
            last_locus = variant_record.chrom

    # write the last record
    if chromosome is not None:
        loc_seqs[chromosome] = record

    # add all remaining loci without variants to the dictionary
    # so that they can be compared with the ground truth
    far = dp.FastaReader(args.stacks_fa)
    for seq, name, *_ in far.chromosomes():

        # Split off the second part of stacks locus names.
        # In the vcf files, the information is not included
        chromosome = name.decode().split()[0]

        # add a record without variants (variant record) for loci without
        # variants detected by stacks.
        if chromosome not in loc_seqs:
            loc_seqs[chromosome] = VCFRecord(chromosome, seq, [], False)

    return list(loc_seqs.values())
