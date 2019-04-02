
rule barcodes:
    output:
        "barcodes/{unit}.tsv"
    run:
        d = individuals.loc[individuals.unit == wildcards.unit,
                            ["p5_barcode", "id"]]
        d[["p5_barcode", "id"]].to_csv(output[0],
                                       index=False, header=None, sep="\t")


rule trim_p7_spacer:
    input:
        lambda w: units.loc[w.unit, "fq2"]
    output:
        "trimmed-spacer/{unit}.2.fq.gz"
    params:
        spacer=lambda w: units.loc[w.unit, "p7_spacer"]
    conda:
        "../envs/seqtk.yaml"
    shell:
        # for b=0, seqtk trimfq uses default behaviour and quality trims
        # to prevent this, only copy the input file, if the spacer length is 0
        """
        if [ {params.spacer} = 0 ]; then
            cp {input} {output}
        else
            seqtk trimfq -b {params.spacer} {input} | gzip > {output}
        fi
        """


rule generate_consensus_reads:
    input:
        fq1=lambda w: units.loc[w.unit, "fq1"],
        fq2="trimmed-spacer/{unit}.2.fq.gz",
    output:
        fq1="dedup/{unit}.consensus.1.fq.gz",
        fq2="dedup/{unit}.consensus.2.fq.gz",
    params:
        umi=config["umi"]
    conda:
        "../envs/consensus.yaml"
    log:
        "logs/consensus/{unit}.log"
    shell:
        "rbt call-consensus-reads -l {params.umi[len]} "
        "-d {params.umi[max_dist]} -D {params.umi[max_seq_dist]} "
        "{input.fq1} {input.fq2} {output.fq1} {output.fq2} 2> {log}"


# remove umi from the p7 read after consensus reads have been computed
rule trim_umi:
    input:
        "dedup/{unit}.consensus.2.fq.gz"
    output:
        "trimmed-umi/{unit}.consensus.2.fq.gz"
    conda:
        "../envs/cutadapt.yaml"
    params:
        umi=config["umi"],
        trim=config["umi"]["len"] + config["restriction-enzyme"]["p7"]["residue-len"]
    log:
        "logs/trim_umi/{unit}.log"
    shell:
        "cutadapt -u {params.trim} {input} -o {output} > {log}"


rule merge_pe_reads:
    input:
        fq1="dedup/{unit}.consensus.1.fq.gz",
        fq2="trimmed-umi/{unit}.consensus.2.fq.gz",
    output:
        merged="merged/{unit}.fq.gz"
    conda:
        "../envs/merge.yaml"
    log:
        "logs/merge/{unit}.log"
    params:
        join_quality='H',
        join_seq=config["reads"]["join_seq"]
    script:
        "../scripts/merge_mates.py"


rule extract:
    input:
        fq1=expand("merged/{unit}.fq.gz", unit=units.id),
        barcodes=expand("barcodes/{unit}.tsv", unit=units.id)
    output:
        expand("extracted/{individual}.fq.gz", individual=individuals.id)
    log:
        expand("logs/extract/{unit}.log",
               unit=units.id)
    params:
        enzymes=config["restriction-enzyme"],
        outdir=get_outdir,
        units=units,
        extra=config["params"]["process_radtags"]
    conda:
        "../envs/stacks.yaml"
    script:
        "../scripts/extract-individuals.py"


rule force_same_length:
    input:
        "extracted/{individual}.fq.gz"
    output:
        "trimmed/{individual}/{individual}.fq.gz"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "len=`seqtk fqchk {input} | grep -oP 'min_len: \\K[0-9]+'`; "
        "seqtk trimfq -L$len {input} | gzip -c > {output}"


rule population_map:
    output:
        "population-map.tsv"
    run:
        d = individuals[["id"]]
        d["pop"] = 1
        d.to_csv(output[0], index=False, header=None, sep="\t")
