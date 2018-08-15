import pysam

header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 10000, 'SN': 'fake'}]}

inbam = pysam.AlignmentFile(snakemake.input[0], "rb", check_sq=False)
outbam = pysam.AlignmentFile(snakemake.output[0], "wb", header=header)

for aln in inbam:
    fake = pysam.AlignedSegment(header=outbam.header)
    m = len(aln.query_sequence)
    # fake.flag = aln.flag | 3
    fake.is_paired = True
    fake.is_proper_pair = True
    fake.is_read1 = aln.is_read1
    fake.is_read2 = aln.is_read2
    fake.cigar = [(0, m)]
    #fake.set_tags(aln.get_tags())
    fake.set_tag("RX", aln.get_tag("RX"))
    fake.set_tag("MQ", 30)
    fake.set_tag("MC", "{}M".format(m))
    fake.mapq = 30
    fake.reference_id = 0
    fake.reference_start = 0
    fake.next_reference_id = 0
    fake.next_reference_start = 0
    fake.query_name = aln.query_name
    fake.query_sequence = aln.query_sequence
    fake.query_qualities = aln.query_qualities
    outbam.write(fake)

inbam.close()
outbam.close()
