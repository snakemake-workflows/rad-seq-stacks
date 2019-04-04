import vcfpy
import yaml
import io

from collections import OrderedDict,namedtuple


SNPRecord = namedtuple("SNPRecord", ["orientation", "pos", "ref", "alt"])
IndelRecord = namedtuple("IndelRecord", ["orientation", "pos", "ref", "alt"])
Allele = namedtuple("Allele", ["cov", "mutations"])

def get_header(individuals):
    header = [
        '##fileformat=VCFv4.3',
        '##fileDate={}'.format("TODAY"),
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        '#' + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + individuals),
        ]
    return io.StringIO("\n".join(header))


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
        # return ("SNP", (orientation, snp_pos, (base_from, base_to)))
        return SNPRecord(orientation, snp_pos, base_from, base_to)
    else:
        # WARNING not implemented yet
        return IndelRecord(None, None, None, None)


def normalize_alleles(individual_alleles, offset):
    normalized_alleles = dict()
    for name, allele in individual_alleles.items():
        normalized_alleles[name] = Allele(
            allele["cov"],
            [normalize_mutation(m, offset) for m in allele["mutations"]]
        )
    return normalized_alleles


def allele_present(individual_alleles, mut):
    """Return a vcf genotype entry:

        '.'   if the individual has a dropout at this locus
        '0|0' if the individual is homozygously ref at the locus
        '1|1' if the individual is either homozygously mutated or
              heterozygously mutated on both alleles
        '0|1' if the individual is heterozygously mutated on one allele

    NOTE:
        This does not work with different mutations of the same position yet.
    """
    if len(individual_alleles) == 0:
        # dropout
        return "."

    elif len(individual_alleles) == 1:
        if 0 in individual_alleles:
            # homozygous reference variant
            return "0|0"
        else:
            # homozygous mutated variant
            return "1|1"

    elif len(individual_alleles) == 2:

        if 0 in individual_alleles:
            # heterozygous mutated variant with reference
            return "0|1"
        else:
            # find out on which allele this variant is present
            on_first, on_second = (0, 0)
            for nr, allele in enumerate(individual_alleles.values()):
                for m in allele.mutations:
                    if mut.pos == m.pos:
                        # a mutation of the active position was found,
                        # now check, whether on the first or on the second
                        # allele
                        if nr == 0:
                            on_first = 1
                        elif nr == 1:
                            on_second = 1
                        else:
                            # simulated genomes are diploid
                            raise ValueError()
            return f"{on_first}|{on_second}"

    else:
        print("This should never happen")
        print(individual_alleles)
        raise ValueError


def generate_records(locus, individuals, chrom):
    records = []
    offset = 100  # KLUDGE

    for allele in locus["allele coverages"].keys():
        if allele == 0:
            continue
        locus_calls = []
        info = OrderedDict()

        # collect all mutated positions
        all_mutations = set()
        for ind in individuals:
            for allele in locus["individuals"][ind].values():
                all_mutations.update(allele["mutations"])

        # normalize them and sort them by position in merged read
        # this is rad seq stacks specific
        normalized_mutations = sorted(
            [
                normalize_mutation(a, offset)
                for a in all_mutations
            ],
            key=lambda mut: mut.pos,
        )

    # TODO: make sure that there is no position with two different alt bases
    # right now, these are not handled properly
    #
    # create one record for each mutation,
    # i.e. each variant at each mutated position
    for mut in normalized_mutations:
        
        # per mutated pos -> record
        # round allele frequencies
        info["AF"] = [round(val, 3) for val in locus["allele frequencies"].values()]
        # coverage of the variant site is the sum of all reads
        # of all individuals
        info["DP"] = sum(locus["allele coverages"].values())

        # check for each individual, if the reference base
        # or another base is present at this location
        for ind in individuals:
            individual_calls = OrderedDict()
            individual_alleles = normalize_alleles(locus["individuals"][ind], offset)
            # coverage for the individual is ths sum of all reads coveraging
            # the site
            individual_calls["DP"] = sum(
                (i.cov for i in individual_alleles.values())
            )
            # get call strings
            individual_calls["GT"] = allele_present(individual_alleles, mut)

            # TODO: handle different variants of the same base on
            # different alleles => REF = A, ALT = C,T, GT= 0|1|2
            locus_calls.append(vcfpy.Call(ind, individual_calls))

        records.append(
            vcfpy.record.Record(
                CHROM=chrom,
                POS=mut.pos,
                ID=[""],
                REF=mut.ref,
                ALT=[vcfpy.Substitution("SNP", mut.alt)],
                QUAL=60,
                FILTER=[],
                INFO=info,
                FORMAT=["GT", "DP"],
                calls=locus_calls
            )
        )
    return records


def parse_rage_gt_file(yaml_path, read_length=100, join_seq="NNNNN",
                       out_path="test.vcf"):
    """TODO
    """
    with open(yaml_path, 'r') as stream:
        try:
            # read all documents in the data
            inds, loci, *other = list(yaml.load_all(stream))
        except yaml.YAMLError as exc:
            print(exc)

    p5_enz = list(inds["Individual Information"].values())[0]["p5 overhang"]
    loc_seqs = []
    # filter out all loci with only one allele, i.e. all unmutated loci
    loci_with_snps = ((n, l) for (n, l) in loci.items()
                      if len(l["allele coverages"]) > 1)

    # print("inds", inds)
    spacer_lengths = [len(i["p5 spacer"]) for i
                      in inds["Individual Information"].values()]
    overhang_lengths = [len(i["p5 overhang"]) for i
                        in inds["Individual Information"].values()]
    offset = read_length - (min(spacer_lengths) + min(overhang_lengths)) \
        + len(join_seq)

    # Assemble header from individual name list
    individuals = [str(i) for i in inds["Individual Information"].keys()]
    reader = vcfpy.Reader.from_stream(
        get_header(individuals))

    with vcfpy.Writer.from_path(out_path, reader.header) as writer:
        for name, locus in loci_with_snps:
            chrom = name.split()[1]

            records = generate_records(locus, individuals, chrom)
            # print("\n".join(map(str, records)))
            for record in records:
                writer.write_record(record)


def main():
    parse_rage_gt_file(yaml_path="../data/easy_dataset/ddRAGEdataset_2_p7_barcodes_gt.yaml")


if __name__ == '__main__':
    main()


# Locus 75:
# allele coverages:
#   0: 86
#   2: 14
#   3: 31
# allele frequencies:
#   0: 0.625
#   2: 0.125
#   3: 0.25
# coverage: 131
# id reads: 0
# individuals:
#   I1:
#     0:
#       cov: 15
#       mutations: [] "0|0"
#     3:
#       cov: 9
#       mutations:
#       - p5@81:G>T  "0|1"
#       - p7@45:G>A  "0|1"
#       - p7@65:C>T  "0|1"
#   I2:
#     0:
#       cov: 38
#       mutations: [] "0|0"
#   I3: {}            "."
#   I4:
#     0:
#       cov: 33
#       mutations: [] "0|0"
#   I5:
#     2:
#       cov: 14
#       mutations:
#       - p7@45:G>A
#       - p7@65:C>T
#     3:
#       cov: 22
#       mutations:
#       - p5@81:G>T "0|1"
#       - p7@45:G>A "1|1"
#       - p7@65:C>T "1|1"
# p5 seq: ATCAGTGCTATGCAGAAGTAGGGTAACAGGCCCTAACATACGCTTACAGATGGTCCTGGTATAGCGCACATGATTATCCCTGGGACGAC
# p7 seq: ACTTGTGACAGTTTAGACCATACGCAAGAGTTAGTGTCCGAATTTGTAACGCCTTCGTGTAATGGCGATCTCCTGATCTAGCAA
