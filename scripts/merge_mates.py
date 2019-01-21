import sys
import click
import dinopy


@click.command(short_help="later")
@click.argument("p5_file", nargs=1, type=click.Path(exists=True))
@click.argument("p7_file", nargs=1, type=click.Path(exists=True))
@click.argument("merged_file", nargs=1, type=click.Path(dir_okay=True))
@click.option("--padding-bases", "padding_bases", default=b"NNNNN")
@click.option("--padding-quality", "padding_quality", default=b"!")
def merge_mates(p5_file, p7_file, merged_file, padding_bases, padding_quality):
    """Merge the p5 and p7 reads from the input files into one read,
    joined by padding bases.
    """
    print(f"Opening files:\n    {p5_file}\n    {p7_file}")
    print(f"Writing to:\n    {merged_file}")
    p5_reader = dinopy.FastqReader(p5_file)
    p7_reader = dinopy.FastqReader(p7_file)

    if len(padding_quality) != 1:
        print("Please specify a single Sanger Phred+33 quality value for "
              "padding.\nExample: --padding-quality !")
        sys.exit(1)
    else:
        padding_quality = padding_quality.encode()

    padding_bases = padding_bases.encode()
    padding_qvs = padding_quality * len(padding_bases)
    with dinopy.FastqWriter(merged_file, force_overwrite=True) as writer:
        for p5_read, p7_read in zip(p5_reader.reads(), p7_reader.reads()):
            merged_seq = b"".join(
                (p5_read.sequence, padding_bases, p7_read.sequence)
            )
            merged_qvs = b"".join(
                (p5_read.quality, padding_qvs, p7_read.quality)
            )
            writer.write(merged_seq, p5_read.name, merged_qvs)


if __name__ == "__main__":
    merge_mates()
