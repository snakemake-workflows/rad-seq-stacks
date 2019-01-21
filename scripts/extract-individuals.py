import os
from snakemake.shell import shell

for barcode_file, fq1, log in zip(snakemake.input.barcodes, snakemake.input.fq1, snakemake.log):
    shell("process_radtags -f {fq1} "
          "--renz_1 {snakemake.params.enzymes[p5]} "
          "-b {barcode_file} -o {snakemake.params.outdir} -y gzfastq -i gzfastq "
          "{snakemake.params.extra} 2> {log}")
    with open(log, "a") as log, open(os.path.join(snakemake.params.outdir, "process_radtags.merged.log")) as stacks_log:
        log.write(stacks_log.read())
