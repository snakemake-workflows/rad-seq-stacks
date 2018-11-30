
rule barcodes:
    output:
        "barcodes/{unit}.tsv"
    run:
        d = individuals.loc[individuals.unit == wildcards.unit, ["p5_barcode", "id"]]
        #d["p7_barcode"] = units.loc[wildcards.unit, "p7_barcode"]
        d[["p5_barcode", "id"]].to_csv(output[0], index=False, header=None, sep="\t")


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


# remove dbr from the p7 read after consensus reads have been computed
rule trim_umi:
    input:
        "dedup/{unit}.consensus.2.fq.gz"
    output:
        "trimmed-umi/{unit}.consensus.2.fq.gz"
    conda:
        "../envs/cutadapt.yaml"
    params:
        dbr="NNNNNNMMGGACG"
    log:
        "logs/trim_umi/{unit}.log"
    shell:
        "cutadapt -g ^{params.dbr} {input} -o {output} 2> {log}"


rule generate_consensus_reads:
    input:
        fq1=lambda w: units.loc[w.unit, "fq1"],
        fq2="trimmed-spacer/{unit}.2.fq.gz",
    output:
        fq1="dedup/{unit}.consensus.1.fq.gz",
        fq2="dedup/{unit}.consensus.2.fq.gz",
    params:
        umi_len=config["params"]["call_consensus_reads"]["umi_len"],
        max_umi_dist=config["params"]["call_consensus_reads"]["max_umi_dist"],
        max_seq_dist=config["params"]["call_consensus_reads"]["max_seq_dist"],
    conda:
        "../envs/consensus.yaml"
    shell:
        "./rbt_release call-consensus-reads -l {params.umi_len} -d {params.max_umi_dist} -D {params.max_seq_dist} {input.fq1} {input.fq2} {output.fq1} {output.fq2}"


rule extract:
    input:
        fq1=expand("dedup/{unit}.consensus.1.fq.gz", unit=units.id),
        fq2=expand("trimmed-umi/{unit}.consensus.2.fq.gz", unit=units.id),
        barcodes=expand("barcodes/{unit}.tsv", unit=units.id)
    output:
        expand(["extracted/{individual}.1.fq.gz",
                "extracted/{individual}.2.fq.gz"],
               individual=individuals.id)
    log:
        expand("logs/extract/{individual}.log",
               individual=individuals.id)
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
        "extracted/{individual}.{read}.fq.gz"
    output:
        "trimmed/{individual}.{read}.fq.gz"
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
