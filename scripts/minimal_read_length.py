"""
"""

lengths = []
for len_file in snakemake.input:
    with open(len_file, 'r') as lf:
        lengths.append(int(lf.readline().strip()))

with open(snakemake.output[0], "w") as out_file:
    print(min(lengths), file=out_file)
