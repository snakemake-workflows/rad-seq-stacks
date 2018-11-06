import csv
import pysam

header = {'HD': {'VN': '1.0', 'SO': 'unsorted', 'GO': 'query', 'SS': 'template-coordinate'}, 'SQ': [{'LN': 10000, 'SN': 'fake'}]}
outbam = pysam.AlignmentFile(snakemake.output[0], "wb", header=header)

def to_bam_rec(cluster_id, dbr, name, seq, qual, read1=True):
    bam_rec = pysam.AlignedSegment(header=outbam.header)
    bam_rec.query_name = name
    bam_rec.query_sequence = seq.encode()
    bam_rec.query_qualities = qual.encode()
    bam_rec.set_tag("RX", dbr.encode())
    bam_rec.set_tag("MI", str(cluster_id))
    bam_rec.cigar = [(0, len(seq))]
    bam_rec.reference_id = 0
    bam_rec.reference_start = 1
    bam_rec.next_reference_id = 0
    bam_rec.next_reference_start = 1
    bam_rec.is_paired = True
    bam_rec.is_proper_pair = True
    bam_rec.is_read1 = read1
    bam_rec.is_read2 = not read1
    outbam.write(bam_rec)

with open(snakemake.input[0]) as clusterfile:
    reader = csv.reader(clusterfile, delimiter="\t")
    for rec in reader:
        cluster_id = rec[0]
        dbr_len = snakemake.config["dbr"]["len"]
        r_seq = rec[7]
        dbr = r_seq[:dbr_len]

        to_bam_rec(cluster_id, dbr, *rec[3:6], read1=True)
        to_bam_rec(cluster_id, dbr, *rec[6:9], read1=False)

outbam.close()
