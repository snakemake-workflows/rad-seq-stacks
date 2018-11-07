import csv
import pysam
import dinopy
import numpy as np
from snakemake.shell import shell


def parse_clusters(stdout):
    for consensus, size, seqids in csv.reader(stdout, delimiter="\t"):
        yield np.fromiter(map(int, seqids.split(",")))


# load dbr sequences
dbrs = np.array([seq[:snakemake.params.dbr_len]
        for _, seq in  dinopy.FastqReader(snakemake.input.fq2)
                             .reads(read_names=False, quality_values=False)])

clusters = dict()

# cluster by read sequences
with subprocess.Popen("starcode --dist {snakemkake.params.seq_dist} --seq-ids "
                      "-1 {snakemake.input.fq1} -2 <(seqtk trimfq "
                      "-b {snakemake.params.dbr_len} {snakemake.input.fq2})",
                      shell=True, stdout=subprocess.PIPE) as seqclust:
    cluster_id = 0
    # iterate over clusters
    for seqids in parse_clusters(seqclust.stdout):
        # get DBRs of clustered sequences
        cluster_dbrs = dbrs[seqids]
        # cluster by DBRs
        with subprocess.Popen("starcode --seq-ids --dist {snakemkake.params.dbr_dist}",
                              shell=True, stdout=subprocess.PIPE) as dbrclust:
            # pass DBRs to cluster process
            for dbr in cluster_dbrs:
                print(dbr, file=dbrclust.stdin)
            # iterate over clusters
            for inner_seqids in parse_clusters(dbrclust.stdout):
                for seqid in seqids[inner_seqids]:
                    clusters[seqid] = cluster_id
                cluster_id += 1


# write clustering results to bam
header = {'HD': {
    'VN': '1.0', 'SO': 'unsorted', 'GO': 'query', 'SS': 'template-coordinate'
    }, 'SQ': [{'LN': 10000, 'SN': 'fake'}]}
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

for seqid, ((f_name, f_seq, f_qual), (r_name, r_seq, r_qual)) in enumerate(zip(
    dinopy.FastqReader(snakemake.input.fq1),
    dinopy.FastqReader(snakemake.input.fq2))):
    cluster_id = clusters[seqid]
    dbr = r_seq[:snakemkake.params.dbr_len]
    to_bam_rec(cluster_id, dbr, f_name, f_seq, f_qual, read1=True)
    to_bam_rec(cluster_id, dbr, r_name, r_seq, r_qual, read2=True)


outbam.close()
