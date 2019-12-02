#!/usr/bin/env python
import subprocess
import os
import glob
import random
import re



def validate_trimmed_spacer():
    """Validate output of trim spacers.

    Invariants:
      - Within one unit, in the p7 reads, all DBRs and all
        enzyme residues should be a the same position within the read.
    """
    print("Evaluating read lengths after spacer trimming (check only first read per unit file):")
    for f in sorted(glob.glob("../trimmed-spacer/*.fq.gz")):
        print()
        for _ in range(10):
            whole_file = subprocess.Popen(["zcat", f], stdout=subprocess.PIPE)
            random_read_start = (random.randint(0, 100000) * 4) + 2
            first_two = subprocess.Popen(('head', f'-n {random_read_start}'), stdin=whole_file.stdout, stdout=subprocess.PIPE)
            output = subprocess.check_output(('tail', '-n 1'), stdin=first_two.stdout)
            print (f"  {f:<30} read in line {random_read_start:>10} length: {len(output):>4}: {output.strip().decode().replace('GGACGTAC', 'GGACG-TAC-', 1)}")


def validate_consensus_reads():
    """Validate output of Generate Consensus Reads.

    Invariants:
      - Within one unit, in the p7 reads, DBRs have been removed, all
        enzyme residues should be a the same position at the beginning
        of the read.
    """
    print("Evaluating the generated consensus reads.")
    p5_enz = re.compile(r"TGCAGG", re.IGNORECASE)
    p7_enz = re.compile(r"^TAC", re.IGNORECASE)

    for f_1, f_2 in zip(sorted(glob.glob("../dedup/*.1.fq.gz")), sorted(glob.glob("../dedup/*.2.fq.gz"))):
        print("p5 -->")
        random_read_starts = [(random.randint(0, 1000) * 4) + 2 for _ in range(10)]
        for start in random_read_starts:
            whole_file = subprocess.Popen(["zcat", f_1], stdout=subprocess.PIPE)
            first_two = subprocess.Popen(('head', f'-n {start}'), stdin=whole_file.stdout, stdout=subprocess.PIPE)
            output = subprocess.check_output(('tail', '-n 1'), stdin=first_two.stdout)
            print (f"  {f_1:<30} p5 read in line {start:>10} length: {len(output):>4}: {re.sub(p5_enz, 'TGCAGG-', output.strip().decode(), 1)}")

        print("<-- p7")
        for start in random_read_starts:
            whole_file = subprocess.Popen(["zcat", f_2], stdout=subprocess.PIPE)
            first_two = subprocess.Popen(('head', f'-n {start}'), stdin=whole_file.stdout, stdout=subprocess.PIPE)
            output = subprocess.check_output(('tail', '-n 1'), stdin=first_two.stdout)
            print (f"  {f_2:<30} p7 read in line {start:>10} length: {len(output):>4}: {re.sub(p7_enz, 'TAC-', output.strip().decode())}")
        print()


def validate_extraction():
    """Validate output of extraction with process_radtags.

    Invariants:
      - More than n percent of the reads are successfully extracted.
        less than 50% -> Error
        less than 80% -> Warning
      - No file should be empty
      - TODO (see below)
    """
    print("Evaluating the extraction of individuals through process_radtags.")

    percentage_line_pattern = re.compile(r"(\d+)(.+)\((.+)%\)")
    def get_percentage(line, pattern):
            result = pattern.search(line.decode())
            return (int(result.group(1)), float(result.group(3)))
        
    for log_file in sorted(glob.glob("../logs/extract/*log")):
        print(f"Analyzing: {log_file}")
        try:
            interesting_lines = subprocess.check_output(["grep", "total sequences", log_file, "-A 4"])
            total, bc_drop, lowq_drop, cutsite_drop, retained = interesting_lines.strip().split(b"\n")
            # bc_total, bc_fraction = get_percentage(bc_drop, percentage_line_pattern)
            # print(f"Dropped {bc_total} reads due to barcodes. {bc_fraction}% of all reads")
            retained_total, retained_fraction = get_percentage(retained, percentage_line_pattern)
            warning_threshold, error_threshold = 95.0, 85.0  # in percent
            if retained_fraction < error_threshold:
                print(f"  WARNING! Very little reads were retained: {retained_fraction}% of reads retained\n")
            elif retained_fraction < warning_threshold:
                print(f"  Some reads were dropped: {retained_fraction}% of reads retained\n")
            else:
                print("  Looking good.\n")
        except subprocess.CalledProcessError as e:
            print(f"  Could not analyze file {log_file}. It might not be finished yet.")
            if e.output == b'':
                print("  No matches found\n")
            else:
                print(e.output + "\n")

    print("Validating file sizes:")
    all_good = True
    for f_1, f_2 in zip(sorted(glob.glob("../extracted/*.1.fq.gz")), sorted(glob.glob("../extracted/*.2.fq.gz"))):
        if "rem" in f_1 or "rem" in f_2:
            # skip files containing removed reads
            continue
        else:
            if os.path.getsize(f_1) == 0:
                print(f"  WARNING: File {f_1} is empty!")
                all_good = False
            if os.path.getsize(f_2) == 0:
                print(f"  WARNING: File {f_2} is empty!")
                all_good = False
    if all_good:
        print("  All good")

    # TODO should read lengths in p5 reads be different after this?
    # TODO are barcodes trimmed out?

def main():
    # validate_trimmed_spacer()
    # validate_consensus_reads()
    validate_extraction()


if __name__ == '__main__':
    main()
