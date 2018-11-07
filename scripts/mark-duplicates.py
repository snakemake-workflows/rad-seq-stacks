import csv
import subprocess

import pysam
import dinopy
import numpy as np

log = open(snakemake.log[0], "w")


def parse_clusters(stdout):
    for consensus, size, seqids in csv.reader(stdout, delimiter="\t"):
        # parse seqids and subtract 1 because starcode provides 1-based indices
        yield np.fromiter(map(int, seqids.split(",")), dtype=int) - 1


# load dbr sequences
dbrs = np.array([seq[:snakemake.params.dbr_len].decode()
        for seq in  dinopy.FastqReader(snakemake.input.fq2)
                          .reads(read_names=False, quality_values=False)])

clusters = dict()

# cluster by read sequences
with subprocess.Popen(f"starcode --dist {snakemake.params.seq_dist} --seq-id "
                      f"-1 <(gzip -d -c {snakemake.input.fq1}) -2 <(seqtk trimfq "
                      f"-b {snakemake.params.dbr_len} {snakemake.input.fq2})",
                      shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                      executable="bash", universal_newlines=True) as seqclust:
    cluster_id = 0
    # iterate over clusters
    for seqids in parse_clusters(seqclust.stdout):
        # get DBRs of clustered sequences
        cluster_dbrs = dbrs[seqids]
        # cluster by DBRs
        with subprocess.Popen(f"starcode --seq-id --dist {snakemake.params.dbr_dist}",
                              shell=True, stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True) as dbrclust:
            # pass DBRs to cluster process
            for dbr in cluster_dbrs:
                print(dbr, file=dbrclust.stdin)
            # close STDIN such that starcode starts to run
            dbrclust.stdin.close()
            # iterate over clusters
            for inner_seqids in parse_clusters(dbrclust.stdout):
                for seqid in seqids[inner_seqids]:
                    clusters[seqid] = cluster_id
                cluster_id += 1
                if (cluster_id - 1) % 1000 == 0:
                    print(f"Processed {cluster_id + 1} clusters.", file=log)
            if dbrclust.wait() != 0:
                raise RuntimeError("Error running starcode: " + dbrclust.stderr)

    if seqclust.wait() != 0:
        raise RuntimeError("Error running starcode: " + seqclust.stderr)

print(f"Total number of clusters: {cluster_id + 1}")


# Write clustering results to bam.
# fgbio requires a sorted BAM file for calculating consensus reads.
# We fake this in the following by letting all reads start at the same position.
header = {'HD': {
    'VN': '1.0', 'SO': 'unsorted', 'GO': 'query', 'SS': 'template-coordinate'
    }, 'SQ': [{'LN': 10000, 'SN': 'fake'}]}
outbam = pysam.AlignmentFile(snakemake.output[0], "wb", header=header)

def to_bam_rec(cluster_id, dbr, fastq_rec, read1=True):
    bam_rec = pysam.AlignedSegment(header=outbam.header)
    bam_rec.query_name = fastq_rec.name
    bam_rec.query_sequence = fastq_rec.sequence
    bam_rec.query_qualities = fastq_rec.quality
    bam_rec.set_tag("RX", dbr)
    bam_rec.set_tag("MI", str(cluster_id))
    bam_rec.cigar = [(0, len(fastq_rec.sequence))]
    bam_rec.reference_id = 0
    bam_rec.reference_start = 0
    bam_rec.next_reference_id = 0
    bam_rec.next_reference_start = 0
    bam_rec.is_paired = True
    bam_rec.is_proper_pair = True
    bam_rec.is_read1 = read1
    bam_rec.is_read2 = not read1
    outbam.write(bam_rec)


# walk over FASTQ records, lookup cluster id, and write everything to a BAM record.
for seqid, (f_rec, r_rec) in enumerate(zip(
    dinopy.FastqReader(snakemake.input.fq1).reads(),
    dinopy.FastqReader(snakemake.input.fq2).reads())):
    cluster_id = clusters[seqid]
    dbr = r_rec.sequence[:snakemake.params.dbr_len]
    to_bam_rec(cluster_id, dbr, f_rec, read1=True)
    to_bam_rec(cluster_id, dbr, r_rec, read1=False)


outbam.close()
log.close()
