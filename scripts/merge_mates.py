"""Merge the p5 and p7 reads from the input files into one read,
joined by join bases.
"""
import sys
import click
import dinopy

# redirect stderr to logfile
sys.stderr = open(snakemake.log[0], "w")

# get and open input files from the calling snakemake rules input directive
p5_file, p7_file = snakemake.input
print(f"Opening files:\n    {p5_file}\n    {p7_file}", file=sys.stderr)
print(f"Writing to:\n    {snakemake.output.merged}", file=sys.stderr)
p5_reader = dinopy.FastqReader(p5_file)
p7_reader = dinopy.FastqReader(p7_file)

# check if the quality value for the join sequence is valid
join_quality = snakemake.params.join_quality
if len(join_quality) != 1:
    print("Please specify a single Sanger Phred+33 quality value for "
          "join_quality.", file=sys.stderr)
    sys.exit(1)
else:
    join_quality = join_quality.encode()

# Assemble join sequence and quality values
join_seq = snakemake.params.join_seq.encode()
join_qvs = (snakemake.params.join_quality * len(join_seq)).encode()

# join all reads
with dinopy.FastqWriter(snakemake.output.merged, force_overwrite=True) as writer:
    for p5_read, p7_read in zip(p5_reader.reads(), p7_reader.reads()):
        merged_seq = b"".join(
            (p5_read.sequence, join_seq, p7_read.sequence)
        )
        merged_qvs = b"".join(
            (p5_read.quality, join_qvs, p7_read.quality)
        )
        writer.write(merged_seq, p5_read.name, merged_qvs)
