"""Take a file simulated with ddrage --multiple-p7-barcodes and split it up by the p7 barcode.

Writes two files named 'reads_<BARCODE_1.fq.gz>' and 'reads_<BARCODE_2.fq.gz>' for each barcode found.
"""
import dinopy as dp
import click


@click.command("plot_single", short_help="Split files simulated by ddRAGE by p7 barcode.")
@click.argument("p5_file", nargs=1, type=click.File(mode='rb'))
@click.argument("p7_file", nargs=1, type=click.File(mode='rb'))
@click.option("--force", "-f", default=False, is_flag=True)
def split_files(p5_file, p7_file, force):
    fqr_fw = dp.FastqReader(p5_file)
    fqr_rev = dp.FastqReader(p7_file)
    output_files = {}
    
    for (fw, rev) in zip(fqr_fw.reads(), fqr_rev.reads()):
        # get bthe nameline of the read (forward or reverse doesn't matter)
        nl = rev.name
        items = nl.split()
    
        # This uses the perfect p7 barcode from the annotation
        # To use the simulated barcode, which can contain sequencing errors,
        # use items[1].split[b":"][-1]
        if items[5].startswith(b"p7_bc"):
            # extract the barcode sequence
            p7_bc = items[5].split(b":")[1].strip(b"'")
    
        # check if a file writer for the barcode is already available
        if p7_bc not in output_files:
            filename_fw = f"reads_{p7_bc.decode()}_1.fq.gz"
            filename_rev = f"reads_{p7_bc.decode()}_2.fq.gz"
            fqw_fw = dp.FastqWriter(filename_fw, force_overwrite=force)
            fqw_rev = dp.FastqWriter(filename_rev, force_overwrite=force)
            fqw_fw.open()
            fqw_rev.open()
            print(f"\nFound new barcode: {p7_bc.decode()}")
            print(f"Writing to:")
            print(f"  -> {filename_fw}")
            print(f"  -> {filename_rev}")
            output_files[p7_bc] = (fqw_fw, fqw_rev)
        else:
            fqw_fw, fqw_rev = output_files[p7_bc]
    
        # write reads back to the writer with the chosen barcode
        fqw_fw.write(*fw)
        fqw_rev.write(*rev)

if __name__ == "__main__":
    split_files()
