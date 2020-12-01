#!/usr/bin/env python
import subprocess
import os
import sys
import random
import re


def get_nr_records_for_fqgz_file(fq_gz_file):
    """Compute number of records in the file using zcat, head and tail.
    """
    whole_file = subprocess.Popen(["zcat", fq_gz_file], stdout=subprocess.PIPE)
    wc_output = subprocess.check_output(
        ('wc', '-l'),
        stdin=whole_file.stdout,
    )
    lines = int(wc_output.strip())
    records = lines // 4
    return records


def extract_read_at(fq_gz_file, pos):
    whole_file = subprocess.Popen(["zcat", fq_gz_file], stdout=subprocess.PIPE)
    random_read_start = (pos * 4) + 2
    first_two = subprocess.Popen(
        ('head', f'-n {random_read_start}'),
        stdin=whole_file.stdout,
        stdout=subprocess.PIPE,
    )
    return subprocess.check_output(
        ('tail', '-n 1'),
        stdin=first_two.stdout,
    )


def validate_trimmed_spacer(out_file, samples_per_file):
    """Validate output of trim spacers.

    NOTE: This test assumes that a DBR is present and is followed
    by the p7 enzyme residue.

    Invariants:
      - Within one unit, in the p7 reads, all DBRs and all
        enzyme residues should be a the same position within the read.
    """
    output_text = [
        "="*120,
        f"Evaluating read lengths after spacer trimming ({samples_per_file} randomly chosen reads per unit file)\n",
        f"Apart from sequencing errors, all of these should show the p7 restriction enzyme in the same position.",
        "="*120,
        ]
    for f in snakemake.input.trimmed_spacer:
        output_text.append(f"Evaluating file: {f:<30}")
        output_text.append(f"{'read in line':>12}\t{'read length':>11}")

        formatted_reads = []
        for _ in range(samples_per_file):

            # randomly pick a record
            records = get_nr_records_for_fqgz_file(f)
            read_nr = random.randint(0, records)
            output = extract_read_at(f, read_nr)

            # Highlight record
            read_seq = output.strip().decode()
            dbr_suffix = snakemake.params.dbr_suffix
            residue = snakemake.params.p7_residue
            highlighted_residue = f"{dbr_suffix}-{residue}-"
            if dbr_suffix + residue in read_seq:
                highlighted_read = read_seq.replace(
                    dbr_suffix + residue,
                    highlighted_residue,
                    1,
                )
            elif dbr_suffix in read_seq:
                highlighted_read = read_seq.replace(
                    dbr_suffix,
                    dbr_suffix+"-",
                    1,
                ) + "-"
            else:
                highlighted_read = read_seq + "--"
            formatted_reads.append(f"{read_nr * 4:>12}\t{len(output):>11} {highlighted_read}")
        # sort reads by position in FASTQ file and append them
        output_text.extend(sorted(formatted_reads))
    output_text.append("\n")
    print("\n".join(output_text), file=out_file)


def validate_consensus_reads(out_file, samples_per_file):
    """Validate output of Generate Consensus Reads.

    Invariants:
      - Within one unit, in the p7 reads, DBRs have been removed, all
        enzyme residues should be a the same position at the beginning
        of the read.
    """
    output_text = [
        "="*120,
        "Evaluating the generated consensus reads. ({samples_per_file} randomly chosen reads per unit file)\n",
        "P5 reads should start with p5 barcode followed by p5 enzyme residue.",
        "P7 reads should no longer contain UMIs/DBRs and begin with the p7 enzyme residue.",
        "="*120,
        ]

    p5_residue = snakemake.params.p5_residue
    p7_residue = snakemake.params.p7_residue
    # Assemble patters to replace restriction enzymes
    p5_enz = re.compile(p5_residue, re.IGNORECASE)
    p7_enz = re.compile(f"^{p7_residue}", re.IGNORECASE)

    for fq_1, fq_2 in zip(snakemake.input.consensus_fq1, snakemake.input.consensus_fq2):
        records = get_nr_records_for_fqgz_file(fq_1)
        random_reads = sorted([random.randint(0, records) for _ in range(samples_per_file)])
        read_nr = random.randint(0, records)

        output_text.append(f"p5 -->   {fq_1:<30}")
        output_text.append(f"{'read in line':>12}\t{'read length':>11}")
        for read_nr in random_reads:
            read_seq = extract_read_at(fq_1, read_nr).strip().decode()
            highlighted_read_seq = re.sub(
                p5_enz,
                f"-{p5_residue}-",
                read_seq,
                1,
            )
            output_text.append(f"{read_nr * 4:>12}\t{len(read_seq):>11}\t{highlighted_read_seq}")

        output_text.append(f"<-- p7   {fq_2:<30}")
        output_text.append(f"{'read in line':>12}\t{'read length':>11}")
        for read_nr in random_reads:
            read_seq = extract_read_at(fq_2, read_nr).strip().decode()
            highlighted_read_seq = re.sub(
                p7_enz,
                p7_residue+"-",
                read_seq,
            )
            output_text.append(f"{read_nr * 4:>12}\t{len(read_seq):>11}\t{highlighted_read_seq}")

        output_text.append("")
    output_text.append("\n")
    print("\n".join(output_text), file=out_file)


def validate_extraction(out_file):
    """Validate output of extraction with process_radtags.

    Invariants:
      - More than n percent of the reads are successfully extracted.
        less than error_threshold -> Error
        less than warning_threshold -> Warning
      - No file should be empty
    """
    output_text = [
        "="*120,
        "Evaluating the extraction of individuals through process_radtags.",
        "="*120,
        ]

    percentage_line_pattern = re.compile(r"(\d+)(.+)\((.+)%\)")

    # Convert to percent points for comparison.
    warning_threshold = snakemake.params.warning_threshold
    error_threshold = snakemake.params.error_threshold

    def get_percentage(line, pattern):
        result = pattern.search(line.decode())
        return (int(result.group(1)), float(result.group(3)))

    raise_error = False
    for log_file in snakemake.input.extract_logs:
        output_text.append(f"Analyzing: {log_file}")

        try:
            interesting_lines = subprocess.check_output(
                ["grep", "total sequences", log_file, "-A 4"]
            )
            total, bc_drop, lowq_drop, cutsite_drop, retained = interesting_lines.strip().split(b"\n")
            retained_total, retained_fraction = get_percentage(
                retained,
                percentage_line_pattern,
            )

            if retained_fraction < error_threshold:
                output_text.append(f"  ERROR! Very little reads were retained: {retained_fraction}% of reads retained\n")
                raise_error = True
            elif retained_fraction < warning_threshold:
                output_text.append(f"  WARNING! Some reads were dropped: {retained_fraction}% of reads retained\n")
            else:
                output_text.append(f"  Looking good: {retained_fraction}% of reads retained\n")
        except subprocess.CalledProcessError as e:
            output_text.append(f"  Could not analyze file {log_file}. It might not be finished yet.")
            if e.output == b'':
                output_text.append("  No matches found\n")
            else:
                output_text.append(e.output + "\n")

    output_text.append("Validating file sizes (check if the file for any individual is empty):")
    all_good = True
    for f_1, f_2 in zip(
            snakemake.input.extracted_fq1,
            snakemake.input.extracted_fq2,
    ):
        if os.path.getsize(f_1) == 0:
            output_text.append(f"  ERROR: File {f_1} is empty!")
            all_good = False
        if os.path.getsize(f_2) == 0:
            output_text.append(f"  ERROR: File {f_2} is empty!")
            all_good = False
    if all_good:
        output_text.append("  All good!")
    else:
        raise_error = True

    if raise_error and snakemake.params.fail_on_error:
        print("\n".join(output_text), file=sys.stderr)
        raise(ValueError("One or more quality thresholds of read_extraction were violated."))
    elif raise_error:
        print(
            "One or more quality thresholds of read_extraction were violated.",
            file=sys.stderr,
        )

    print("\n".join(output_text), file=out_file)


def main():
    samples_per_file = snakemake.params.samples_per_file
    import sys

    with open(snakemake.log[0], "w") as log_file:
        # redirect all output top the logfile as suggested
        # here: https://stackoverflow.com/a/64114743/2862719
        sys.stderr = log_file
        with open(snakemake.output.qc_log, "w") as out_file:
            validate_trimmed_spacer(out_file, samples_per_file)
            validate_consensus_reads(out_file, samples_per_file)
            validate_extraction(out_file)


if __name__ == '__main__':
    main()
