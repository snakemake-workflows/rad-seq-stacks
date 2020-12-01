"""Translate a ddRAGE ground truth YAML file into vcf format.

This workflow is adapted to the rad-seq-stacks pipeline and
uses assumptions specific to this pipeline, most notably
that p5 and p7 reads are concatenated before analysis.
"""
import vcfpy
import yaml
import io
import datetime
import sys

from collections import OrderedDict, namedtuple, defaultdict


SNPRecord = namedtuple("SNPRecord", ["orientation", "pos", "ref", "alt"])
IndelRecord = namedtuple("IndelRecord", ["orientation", "pos", "ref", "alt"])
Allele = namedtuple("Allele", ["cov", "mutations"])


def get_header(individuals):
    """Format a VCF header, including individual names and the current date."""
    time = datetime.datetime.now()
    header = [
        '##fileformat=VCFv4.3',
        f'##fileDate={time.year}{time.month:0>2}{time.day:0>2}',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples '
        'With Data">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele Depth">',
        '##FORMAT=<ID=GL,Number=G,Type=Float,Description='
        '"Genotype Likelihood">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description='
        '"Genotype Quality">',
        '#' + "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                         "INFO", "FORMAT"] + individuals),
        ]
    return io.StringIO("\n".join(header))


def parse_mutation(mut, offset):
    """Parse a ddRAGE mutation into a SNPRecord object.

    If a mutation falls on the p7 read, add the offset so that this works
    with concatenated reads.

    Returns:
        NamedTuple (SNPRecord or IndelRecord)
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
        return SNPRecord(orientation, snp_pos, base_from, base_to)
    else:
        # WARNING not implemented yet
        return IndelRecord(None, None, None, None)


def parse_alleles(individual_alleles, offset):
    """Parse all mutations of the alleles."""
    normalized_alleles = dict()
    for name, allele in individual_alleles.items():
        normalized_alleles[name] = Allele(
            allele["cov"],
            [parse_mutation(m, offset) for m in allele["mutations"]]
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
        return None, "."

    elif len(individual_alleles) == 1:
        if 0 in individual_alleles:
            # homozygous reference variant
            return (0, 0), "0|0"
        else:
            # homozygous mutated variant
            # check if the homozygous variant is the one in question
            for allele in individual_alleles.values():
                for m in allele.mutations:
                    if mut.pos == m.pos:
                        return (1, 1), "1|1"
            # if this is reached, the mutation in question was not present
            return (0, 0), "0|0"

    elif len(individual_alleles) == 2:

        if 0 in individual_alleles:
            # heterozygous mutated variant with reference
            return (0, 1), "0|1"
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
            return (on_first, on_second), f"{on_first}|{on_second}"

    else:
        print("Invalid number of individual alleles. "
              "No individual should have more than two alleles!",
              file=sys.stderr)
        print("Alleles dictionary:", individual_alleles, file=sys.stderr)
        raise ValueError


def generate_records(locus, individuals, chrom, offset):
    """Generate VCF records for all mutations of the given locus."""
    records = []

    mutation_allele_frequencies = dict()
    mutation_sample_count = defaultdict(set)

    # assemble normalized and sorted list of mutations, each of which
    # will later receive one record in the VCF file
    # additionally, create a dictionary that maps mutations to their
    # allele frequency: the sum of all alleles with this specific mutation
    for allele, frequency in locus["allele frequencies"].items():
        # print("allele", allele)
        if allele == 0:
            # skip reference alleles
            continue
        else:
            # collect all mutated positions
            all_mutations = set()
            for ind in individuals:
                mutations_per_individual = set()
                try:
                    mutations_per_individual.update(
                        locus["individuals"][ind][allele]["mutations"])
                    # print(mutations_per_individual)
                    for m in mutations_per_individual:
                        mutation_sample_count[parse_mutation(m, offset)].add
                        (ind)
                    all_mutations.update(mutations_per_individual)
                except KeyError:
                    # this allele is not present in this individual
                    ...

            # normalize them and sort them by position in merged read
            # this is rad seq stacks specific
            normalized_mutations = sorted(
                [
                    parse_mutation(a, offset)
                    for a in all_mutations
                ],
                key=lambda mut: mut.pos,
            )

            for mut in normalized_mutations:
                if mut in mutation_allele_frequencies:
                    mutation_allele_frequencies[mut] += frequency
                else:
                    mutation_allele_frequencies[mut] = frequency

    individual_allele_coverage = defaultdict(lambda: 0)
    for ind_name, alleles in locus["individuals"].items():
        for name, allele in alleles.items():
            if not allele["mutations"]:
                individual_allele_coverage[ind_name, 0] += allele["cov"]
            else:
                for mut in allele["mutations"]:
                    individual_allele_coverage[ind_name, parse_mutation(mut, offset)] += allele["cov"]

    # TODO: make sure that there is no position with two different alt bases
    # right now, these are not handled properly
    #
    # create one record for each mutation,
    # i.e. each variant at each mutated position
    # print("norm mut", normalized_mutations)
    for mut in normalized_mutations:
        info = OrderedDict()
        locus_calls = []
        # print(f"looking for {mut}")
        # per mutated pos -> record

        # round allele frequencies
        info["AF"] = [round(mutation_allele_frequencies[mut], 3)]
        # coverage of the variant site is the sum of all reads
        # of all individuals
        info["DP"] = sum(locus["allele coverages"].values())
        # number of samples with mutation
        info["NS"] = len(mutation_sample_count[mut])

        # check for each individual, if the reference base
        # or another base is present at this location
        for ind in individuals:
            individual_calls = OrderedDict()
            individual_alleles = parse_alleles(locus["individuals"][ind],
                                               offset)
            # print(f"norm individual alleles for {ind}", individual_alleles)
            # coverage for the individual is ths sum of all reads coveraging
            # the site
            individual_calls["DP"] = sum(
                (i.cov for i in individual_alleles.values())
            )
            # get call strings
            allele_presence, allele_str = allele_present(individual_alleles, mut)
            individual_calls["GT"] = allele_str
            # print(individual_calls["GT"])

            # fill individual allele coverage as a tuple of
            # (coverage of ref allele, coverage of alt allele)
            try:
                
                if allele_presence is None:
                    ind_allele_cov = (0, 0)
                elif allele_presence == (0, 0):
                    ind_allele_cov = (individual_allele_coverage[(ind, 0)], 0)
                elif allele_presence == (1, 1):
                    ind_allele_cov = (0, individual_allele_coverage[(ind, mut)])
                elif allele_presence == (1, 0):
                    ind_allele_cov = (individual_allele_coverage[(ind, mut)], individual_allele_coverage[(ind, 0)])
                elif allele_presence == (0, 1):
                    ind_allele_cov = (individual_allele_coverage[(ind, 0)], individual_allele_coverage[(ind, mut)])
                else:
                    raise ValueError("Invalid mutation")
                individual_calls["AD"] = ind_allele_cov
            except KeyError:
                print(individual_allele_coverage)
                raise
            
            # TODO: handle different variants of the same base on
            # different alleles => REF = A, ALT = C,T, GT= 0|1|2
            locus_calls.append(vcfpy.Call(ind, individual_calls))

        rec = vcfpy.record.Record(
                CHROM=chrom,
                POS=mut.pos,
                ID=[""],
                REF=mut.ref,
                ALT=[vcfpy.Substitution("SNP", mut.alt)],
                QUAL="",
                FILTER=["PASS"],
                INFO=info,
                FORMAT=["GT", "DP", "AD"],
                calls=locus_calls
            )
        # print("Record:", rec)
        records.append(
            rec
        )
    return records


def parse_rage_gt_file(yaml_path, read_length, join_seq, out_path):
    """Create a VCF file that contains one record for each mutation simulated
    by ddRAGE. Loci with multiple simulated mutations receive multiple lines.
    """
    # Read ground truth file. Skip HRL reads and singletons after
    # the loci, if those are present.
    with open(yaml_path, 'r') as stream:
        try:
            # read all documents in the data
            inds, loci, *other = list(yaml.load_all(stream, Loader=yaml.FullLoader))
        except yaml.YAMLError as exc:
            print(exc)

    # filter out all loci with only one allele, i.e. all unmutated loci
    loci_with_snps = ((n, l) for (n, l) in loci.items()
                      if len(l["allele coverages"]) > 1)

    # Compute the offset for mutations on the concatenated p7 read.
    # The added length is the total simulated read size minus the shortest
    # possible adapter sequence (i.e. the longest possible genomic sequence)
    # plus the length of the sequence of 'N's used to join the p5 and p7 reads.
    spacer_lengths = [len(i["p5 spacer"]) for i
                      in inds["Individual Information"].values()]
    overhang_lengths = [len(i["p5 overhang"]) for i
                        in inds["Individual Information"].values()]
    offset = read_length - (min(spacer_lengths) + min(overhang_lengths)) \
        + len(join_seq)

    # Assemble VCF header from individual name list
    individuals = [str(i) for i in inds["Individual Information"].keys()]
    reader = vcfpy.Reader.from_stream(
        get_header(individuals))

    # Iterate through all ddRAGE loci and create a VCF record for each
    # simulated mutation
    with vcfpy.Writer.from_path(out_path, reader.header) as writer:
        for name, locus in loci_with_snps:
            # Get the locus number, by splitting off the leading 'Locus'
            chrom = name.split()[1]
            # build VCF records and write them to file
            for record in generate_records(locus, individuals, chrom, offset):
                writer.write_record(record)


parse_rage_gt_file(
    yaml_path=snakemake.input.yaml,  # noqa: F821
    read_length=snakemake.params.read_length,  # noqa: F821
    join_seq=snakemake.params.join_seq,  # noqa: F821
    out_path=snakemake.output.vcf,  # noqa: F821
)
