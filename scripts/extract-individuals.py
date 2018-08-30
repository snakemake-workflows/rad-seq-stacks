from snakemake.shell import shell

for barcode_file, fq1, fq2 in zip(snakemake.input.barcodes, snakemake.input.fq1, snakemake.input.fq2):
    shell("process_radtags -1 {fq1} -2 {fq2} "
          "--renz_1 {snakemake.params.enzymes[p5]} "
          "--renz_2 {snakemake.params.enzymes[p7]} "
          "-b {barcode_file} -o {snakemake.params.outdir} -y gzfastq -i gzfastq "
          "{snakemake.params.extra}")
