rule shortest_read_per_sample:
    input:
        "{path}/{individual}.{orientation}.fq.gz"
    output:
        "{path}/{individual}.{orientation}.len"
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk fqchk {input} | grep -oP 'min_len: \\K[0-9]+' > {output}"

rule shortest_read:
    input:
        expand("{{path}}/{individual}.{{orientation}}.len", individual=individuals.id)
    output:
        "{path}/all.{orientation}.len"
    script:
        "../scripts/minimal_read_length.py"

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
    benchmark:
        "benchmarks/trimmed-spacer/{unit}.txt"
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
    benchmark:
        "benchmarks/consensus/{unit}.txt"
    shell:
        "TMPDIR=dedup "
        "rbt call-consensus-reads -l {params.umi[len]} --umi-on-reverse "
        "-d {params.umi[max_dist]} -D {params.umi[max_seq_dist]} "
        "{input.fq1} {input.fq2} {output.fq1} {output.fq2} 2> {log}"


# remove restriction enzyme residue p7 read after individual extraction
# consensus reads have already been computed and umis have already been
# removed during consensus read generation, so the residue is at the start
# of the p7 read
rule trim_residue:
    input:
        "extracted/{individual}.2.fq.gz"
    output:
        "trimmed-residue/{individual}.2.fq.gz"
    conda:
        "../envs/cutadapt.yaml"
    params:
        trim=config["restriction-enzyme"]["p7"]["residue-len"]
    log:
        "logs/trim_residue/{individual}.log"
    benchmark:
        "benchmarks/trimmed-residue/{individual}.txt"
    shell:
        "cutadapt -u {params.trim} {input} -o {output} > {log}"


rule merge_pe_reads:
    input:
        fq1="extracted/{individual}.1.fq.gz",
        fq2="trimmed-residue/{individual}.2.fq.gz",
        fq1_length="extracted/all.1.len",
        fq2_length="trimmed-residue/all.2.len",
    output:
        merged="merged/{individual}.fq.gz"
    conda:
        "../envs/merge.yaml"
    log:
        "logs/merge/{individual}.log"
    benchmark:
        "benchmarks/merge/{individual}.txt"
    params:
        join_quality='H',
        join_seq=config["reads"]["join_seq"]
    script:
        "../scripts/merge_mates.py"


# concatenate p5 and p7 files into one .fq.gz file
rule concatenate_read_files:
    input:
        fq1="extracted/{individual}.1.fq.gz",
        fq2="trimmed-residue/{individual}.2.fq.gz",
    output:
        concat="concatenated/{individual}.fq.gz"
    log:
        "logs/concatenate/{individual}.log"
    benchmark:
        "benchmarks/merge/{individual}.txt"
    shell:
        "zcat {input.fq1} {input.fq2} | gzip -c > {output.concat}"


rule extract:
    input:
        fq1=expand("dedup/{unit}.consensus.1.fq.gz", unit=units.id),
        fq2=expand("dedup/{unit}.consensus.2.fq.gz", unit=units.id),
        barcodes=expand("barcodes/{unit}.tsv", unit=units.id)
    output:
        expand(["extracted/{individual}.1.fq.gz",
                "extracted/{individual}.2.fq.gz"],
               individual=individuals.id)
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



def trim_input():
    """Depending on the desired mode, this can either be p5 reads only,
    p5 and p7 reads merged into one ('single end') read, or the p5 and p7 
    files concatenated into one file.
    """
    mode = config["reads"]["mode"]
    if mode == "p5_only":
        return "extracted/{individual}.1.fq.gz"
    elif mode == "merged":
        return "merged/{individual}.fq.gz"
    elif mode == "concatenated":
        return "concatenated/{individual}.fq.gz"
    else:
        raise ValueError(f"Invalid mode: {mode}. Should be 'p5_only', 'merged', or 'concatenated'")


# Trim all (merged) reads of one individual to the same length
rule force_same_length:
    input:
        trim_input()
    output:
        "trimmed/{individual}/{individual}.fq.gz"
    benchmark:
        "benchmarks/trim_lenth/{individual}.txt"
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
