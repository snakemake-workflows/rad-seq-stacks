import dinopy
import sys

in_fastq = dinopy.FastqReader(sys.stdin)
with dinopy.FastqWriter(sys.stdout) as out_fastq:
    out_fastq.write_reads((read.sequence, read.name.split(b" ")[0], read.quality) 
                          for read in in_fastq.reads())

