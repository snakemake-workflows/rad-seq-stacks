rule shortest_read_per_sample:
    input:
        "{path}/{individual}.{orientation}.fq.gz"
    output:
        "{path}/{individual}.{orientation}.len"
    conda:
        "../envs/seqtk.yaml"
    log:
        "logs/shortest_read_per_sample/{path}/{individual}.{orientation}.log"
    shell:
        "seqtk fqchk {input} 2> {log} | grep -oP 'min_len: \\K[0-9]+' > {output} 2> {log}"


rule shortest_read:
    input:
        expand("{{path}}/{individual}.{{orientation}}.len",
               individual=individuals.id)
    output:
        "{path}/all.{orientation}.len"
    log:
        "logs/shortest_read/{path}/{orientation}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/minimal_read_length.py"


rule barcodes:
    output:
        "resources/barcodes/{unit}.tsv"
    log:
        "logs/barcodes/{unit}.log"
    params:
        individuals=individuals
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/barcodes.py"


rule trim_p7_spacer:
    input:
        lambda w: units.loc[w.unit, "fq2"]
    output:
        "analysis/trimmed-spacer/{unit}.2.fq.gz"
    benchmark:
        "benchmarks/trimmed-spacer/{unit}.txt"
    params:
        spacer=lambda w: units.loc[w.unit, "p7_spacer"]
    conda:
        "../envs/seqtk.yaml"
    log:
        "logs/trim_p7_spacer/{unit}.log"
    shell:
        # for b=0, seqtk trimfq uses default behaviour and quality trims
        # to prevent this, only copy the input file, if the spacer length is 0
        """
        if [ {params.spacer} = 0 ]; then
            cp {input} {output}
        else
            seqtk trimfq -b {params.spacer} {input} 2> {log} | gzip > {output} 2> {log}
        fi
        """


rule generate_consensus_reads:
    input:
        fq1=lambda w: units.loc[w.unit, "fq1"],
        fq2="analysis/trimmed-spacer/{unit}.2.fq.gz",
    output:
        fq1="analysis/dedup/{unit}.consensus.1.fq.gz",
        fq2="analysis/dedup/{unit}.consensus.2.fq.gz",
    params:
        umi=config["umi"],
        tmpdir=lambda w, output: os.path.dirname(output[0])
    conda:
        "../envs/consensus.yaml"
    log:
        "logs/consensus/{unit}.log"
    benchmark:
        "benchmarks/consensus/{unit}.txt"
    shell:
        "TMPDIR={params.tmpdir} "
        "rbt call-consensus-reads -l {params.umi[len]} --umi-on-reverse "
        "-d {params.umi[max_dist]} -D {params.umi[max_seq_dist]} "
        "{input.fq1} {input.fq2} {output.fq1} {output.fq2} 2> {log}"


# remove restriction enzyme residue p7 read after individual extraction
# consensus reads have already been computed and umis have already been
# removed during consensus read generation, so the residue is at the start
# of the p7 read
rule trim_residue:
    input:
        "analysis/extracted/{individual}.2.fq.gz"
    output:
        "analysis/trimmed-residue/{individual}.2.fq.gz"
    conda:
        "../envs/cutadapt.yaml"
    params:
        trim=config["restriction-enzyme"]["p7"]["residue-len"]
    log:
        "logs/trim_residue/{individual}.log"
    benchmark:
        "benchmarks/trimmed-residue/{individual}.txt"
    shell:
        "cutadapt -u {params.trim} {input} -o {output} &> {log}"


rule replace_residue:
    input:
        "analysis/extracted/{individual}.2.fq.gz"
    output:
        "analysis/replaced-residue/{individual}.2.fq.gz"
    log:
        "logs/replace_residue/{individual}.log"
    benchmark:
        "benchmarks/replace-residue/{individual}.txt"
    params:
        p5_res=config["restriction-enzyme"]["p5"]["residue-seq"],
        p7_res=config["restriction-enzyme"]["p7"]["residue-seq"],
    conda:
        "../envs/sed.yaml"
    shell:
        "zcat {input} 2> {log} | "
        "sed -e 's/^{params.p7_res}/{params.p5_res}/' 2> {log} | "
        "gzip -c > {output} 2> {log}"


rule merge_pe_reads:
    input:
        fq1="analysis/extracted/{individual}.1.fq.gz",
        fq2="analysis/trimmed-residue/{individual}.2.fq.gz",
        fq1_length="analysis/extracted/all.1.len",
        fq2_length="analysis/trimmed-residue/all.2.len",
    output:
        merged="analysis/merged/{individual}.fq.gz"
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
        fq1="analysis/extracted/{individual}.1.fq.gz",
        fq2="analysis/replaced-residue/{individual}.2.fq.gz",
    output:
        concat="analysis/concatenated/{individual}.fq.gz"
    log:
        "logs/concatenate/{individual}.log"
    benchmark:
        "benchmarks/merge/{individual}.txt"
    conda:
        "../envs/python.yaml"
    shell:
        "zcat {input.fq1} {input.fq2} 2> {log} | gzip -c > {output.concat} 2> {log}"


rule extract:
    input:
        fq1=expand("analysis/dedup/{unit}.consensus.1.fq.gz", unit=units.unit),
        fq2=expand("analysis/dedup/{unit}.consensus.2.fq.gz", unit=units.unit),
        barcodes=expand("resources/barcodes/{unit}.tsv", unit=units.unit)
    output:
        expand(["analysis/extracted/{individual}.1.fq.gz",
                "analysis/extracted/{individual}.2.fq.gz"],
               individual=individuals.id)
    log:
        expand("logs/extract/{unit}.log",
               unit=units.unit)
    params:
        enzymes=config["restriction-enzyme"],
        outdir=get_outdir,
        units=units,
        extra=config["params"]["process_radtags"]
    conda:
        "../envs/stacks.yaml"
    script:
        "../scripts/extract-individuals.py"


# Trim all (merged) reads of one individual to the same length
rule force_same_length:
    input:
        trim_input()
    output:
        "analysis/trimmed/{individual}/{individual}.fq.gz"
    benchmark:
        "benchmarks/trim_lenth/{individual}.txt"
    conda:
        "../envs/seqtk.yaml"
    log:
        "logs/force_same_length/{individual}.log"
    shell:
        "len=`seqtk fqchk {input} 2> {log} | grep -oP 'min_len: \\K[0-9]+'`; "
        "seqtk trimfq -L$len {input} 2> {log} | gzip -c > {output} 2> {log}"


rule population_map:
    output:
        "resources/population-map.tsv"
    log:
        "logs/population_map/population_map.log"
    params:
        individuals=individuals
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/population_map.py"
