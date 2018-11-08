"""Use ddRAGE ground truth files to validate the loci assembled by stacks.

Maybe this should be moved into the .test folder.
"""

rule index_stacks_loci:
    input:
        "stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.fa.gz"
    output:
        fa="stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.fa",
        fai="stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.fa.fai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        gunzip -k {input}
        samtools faidx {output.fa}
        """


rule compare_loci:
    input:
        fa="stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.fa",
        fai="stacks/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/catalog.fa.fai",
        vcf="calls/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/populations.snps.vcf",
        gt="data/raw_dataset/ddRAGEdataset_2_p7_barcodes_gt.yaml",
    output:
        "validation/n={max_locus_mm}.M={max_individual_mm}.m={min_reads}/validation.txt"
    conda:
        "../envs/testdata.yaml"
    shell:
        "python ../scripts/evaluate_stacks_results.py --ground-truth {input.gt} --stacks-snps-file {input.vcf} --stacks-fasta-file {input.fa} -o {output}"
