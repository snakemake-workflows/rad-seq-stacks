import os
from snakemake.shell import shell

for barcode_file, fq1, fq2, log in zip(snakemake.input.barcodes, snakemake.input.fq1, snakemake.input.fq2, snakemake.log):
    shell("process_radtags -1 {fq1} -2 {fq2} "
          "--renz_1 {snakemake.params.enzymes[p5][name]} "
          "--renz_2 {snakemake.params.enzymes[p7][name]} "
          "-b {barcode_file} -o {snakemake.params.outdir} -y gzfastq -i gzfastq "
          "{snakemake.params.extra} 2> {log}")
    with open(log, "a") as log, open(os.path.join(snakemake.params.outdir, "process_radtags.dedup.log")) as stacks_log:
        log.write(stacks_log.read())
