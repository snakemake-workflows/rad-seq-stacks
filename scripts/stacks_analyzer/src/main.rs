//! Take a stacks .tags.tsv file generated using only p5 reads and assemble the
//! corresponding p7 reads. This is to find out what is going on on the p7 side.
//!
use bio::io::fastq;
use flate2::bufread::MultiGzDecoder;
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::BufWriter;
use structopt::StructOpt;

/// Create a read_name to PE read (fastq records) mapping from a pair of
/// .fastq.gz files
fn get_mapping(
    p5_path: &str,
    p7_path: &str,
) -> HashMap<String, (bio::io::fastq::Record, bio::io::fastq::Record)> {
    // Set up progress spinner
    let spinner_style = indicatif::ProgressStyle::default_spinner()
        .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ")
        .template("{prefix:.bold.dim} {spinner} {wide_msg}");
    let pb = indicatif::ProgressBar::new_spinner();
    pb.set_style(spinner_style.clone());
    pb.set_prefix(&format!(
        "[1/3] Parsing read files {} and {}",
        p5_path, p7_path
    ));

    // Creater readers
    let p5_reader = fastq::Reader::new(
        fs::File::open(p5_path)
            .map(BufReader::new)
            .map(MultiGzDecoder::new)
            .unwrap(),
    );
    let p7_reader = fastq::Reader::new(
        fs::File::open(p7_path)
            .map(BufReader::new)
            .map(MultiGzDecoder::new)
            .unwrap(),
    );

    // Read in pairs of fastq records and store them in a hash map
    let mut name_to_pe = HashMap::new();
    for (i, (p5_read, p7_read)) in p5_reader.records().zip(p7_reader.records()).enumerate() {
        let p5_read = p5_read.unwrap();
        let p7_read = p7_read.unwrap();
        name_to_pe.insert(
            String::from(p5_read.id()),
            (p5_read.clone(), p7_read.clone()),
        );
        // tend to progress spinner
        pb.inc(1);
        if i % 1_000 == 0 {
            pb.set_message(&format!("  Processed {:>10} reads", i));
        };
    }
    pb.finish_with_message(&"Done.".to_string());

    name_to_pe
}

/// Struct to model one entry in a stacks .tags.tsv file
/// Right now contains more info than is actually needed.
#[derive(Debug)]
struct StacksCluster {
    name: String,
    primary_reads: Vec<(String, String)>,
    secondary_reads: Vec<(String, String)>,
    consensus: Option<String>,
    model_line: Option<String>,
}

impl StacksCluster {
    /// Build a stacks cluster from a .tags.tsv file.
    /// Write the length of each parsed cluster to the dat file to plot
    /// them later
    fn from_csv_records(tsv_path: &str, data_path: &str) -> Vec<Self> {
        // set up preogress spinner
        let spinner_style = indicatif::ProgressStyle::default_spinner()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ")
            .template("{prefix:.bold.dim} {spinner} {wide_msg}");
        let pb = indicatif::ProgressBar::new_spinner();
        pb.set_style(spinner_style.clone());
        pb.set_prefix(&format!(
            "[2/3] Parsing stacks .tags.tsv file: {}",
            tsv_path
        ));

        // Prepare input and output files
        let tsv_file = fs::File::open(tsv_path)
            .map(BufReader::new)
            .map(MultiGzDecoder::new)
            .unwrap();
        let mut data_file = match File::create(data_path) {
            Err(e) => panic!("Couldn't write dat file: {:?}", e),
            Ok(data_file) => BufWriter::new(data_file),
        };

        let mut primary = vec![];
        let mut secondary = vec![];
        let mut stacks = vec![];
        let mut name: Option<String> = None;
        let mut model: Option<String> = None;
        let mut consensus: Option<String> = None;
        let mut nr_cluster = 0;

        // Parse entries of .tags.tsv file
        // Each stack comprises:
        // - consensus line containing the stack consensus sequence
        // - model line with mutation information for each position
        // - several primary reads
        // - optionally several secondary reads
        for rec in csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(true)
            .from_reader(tsv_file)
            .records()
        {
            let rec = rec.unwrap();

            // find put which line type is present
            match rec[2].as_ref() {
                "consensus" => {
                    // Consensus lines start a new stack.
                    // Write size data to dat file
                    // Create a new StacksCluster object with the collected
                    // reads and keep it in memory.
                    match consensus {
                        Some(c) => {
                            // Finalize an existing cluster

                            // tend to the progress spinner
                            nr_cluster += 1;
                            pb.inc(1);
                            if nr_cluster % 1_000 == 0 {
                                pb.set_message(&format!("  Processed {:>10} cluster", nr_cluster));
                            };

                            // Write size data to dat file
                            data_file
                                .write_all(
                                    format!("{}\n", primary.len() + secondary.len()).as_bytes(),
                                )
                                .unwrap();

                            // Build stack and keep it in stacks list
                            // empty both read vectors by draining them
                            stacks.push(StacksCluster {
                                name: name.clone().unwrap(),
                                consensus: Some(c),
                                model_line: model.clone(),
                                primary_reads: primary.drain(0..).collect(),
                                secondary_reads: secondary.drain(0..).collect(),
                            });
                            // Set consensus line to the consensus of the new stack
                            consensus = Some(String::from(&rec[1]));
                        }
                        None => {
                            // Initialize consensus line for the first stack
                            consensus = Some(String::from(&rec[1]));
                        }
                    }
                }
                "model" => {
                    // Collect name and model line information for the stack
                    name = Some(String::from(&rec[1]));
                    model = Some(String::from(&rec[5]));
                }
                "primary" => {
                    // Add a primary read name and sequence to the stack
                    primary.push((String::from(&rec[4]), String::from(&rec[5])));
                }
                "secondary" => {
                    // Add a secondary read name and sequence to the stack
                    secondary.push((String::from(&rec[4]), String::from(&rec[5])));
                }
                _ => {
                    // Line couldn't be parsed by one of the previous categories.
                    // This should not happen in a valid .tags.tsv file
                    println!("{:?}", rec);
                    panic!("Could not parse record!");
                }
            }
        }
        // finish last stack
        if let (Some(name), Some(consensus), Some(model)) = (name, consensus, model) {
            data_file
                .write_all(format!("{}\n", primary.len() + secondary.len()).as_bytes())
                .unwrap();

            // Build stack and keep it in stacks list
            // empty both read vectors by draining them
            stacks.push(StacksCluster {
                name,
                consensus: Some(consensus),
                model_line: Some(model),
                primary_reads: primary.drain(0..).collect(),
                secondary_reads: secondary.drain(0..).collect(),
            });
        }

        pb.finish_with_message(&format!("Done. Assembled {:} cluster.", nr_cluster));
        stacks
    }
}

/// Use mapping and stacks cluster to compile p5+p7 pileups
fn assemble_clusters(
    mapping: HashMap<String, (bio::io::fastq::Record, bio::io::fastq::Record)>,
    stacks_cluster: Vec<StacksCluster>,
    out_path: &str,
) {
    // create progress bar
    let bar_style = indicatif::ProgressStyle::default_bar()
        .template(&format!("{{prefix:.bold}}▕{{bar:.{}}}▏{{msg}}", "blue"))
        .progress_chars("█▉▊▋▌▍▎▏  ");
    let length = stacks_cluster.len();
    let pb = indicatif::ProgressBar::new(length as u64);
    pb.set_style(bar_style.clone());
    pb.set_prefix(&"[3/3] Assembling Loci and writing to stdout.".to_string());

    // Prepare output file
    let mut out_file = match File::create(out_path) {
        Err(e) => panic!("Couldn't write output file: {:?}", e),
        Ok(data_file) => BufWriter::new(data_file),
    };

    let mut primary_seqs = Vec::with_capacity(100);
    let mut secondary_seqs = Vec::with_capacity(100);

    // Iterate through all reads in a stack, collect the corresponding p5 p7 pair
    // from the mapping (i.e. from the FASTQ files) and write the pileup to file
    for (i, cluster) in stacks_cluster.iter().enumerate() {
        out_file
            .write_all(format!("\n\nCluster {:}\n\n", i).as_bytes())
            .unwrap();
        // Collect all primary reads into a vector
        for read in cluster.primary_reads.iter() {
            let (name, _seq) = read;
            let (p5_rec, p7_rec) = mapping.get(name).unwrap();
            primary_seqs.push((
                std::str::from_utf8(p5_rec.seq()).unwrap(),
                std::str::from_utf8(p7_rec.seq()).unwrap(),
            ));
        }
        // Collect all secondary reads into a different vector so that they
        // can be handled independently from primary reads
        for read in cluster.secondary_reads.iter() {
            let (name, _seq) = read;
            let (p5_rec, p7_rec) = mapping.get(name).unwrap();
            secondary_seqs.push((
                std::str::from_utf8(p5_rec.seq()).unwrap(),
                std::str::from_utf8(p7_rec.seq()).unwrap(),
            ));
        }
        // Sort read pairs by p7 read. p5 reads should already be similar enough
        // by being in the same stack in the tsv file
        primary_seqs.sort_by(|r1, r2| r1.1.cmp(r2.1));
        secondary_seqs.sort_by(|r1, r2| r1.1.cmp(r2.1));
        // Write primaries to file
        for (p5, p7) in primary_seqs.drain(0..) {
            out_file
                .write_all(format!("{}---{}\n", p5, p7).as_bytes())
                .unwrap();
        }
        out_file.write_all("---\n".to_string().as_bytes()).unwrap();
        // Write secondaries to file
        for (p5, p7) in secondary_seqs.drain(0..) {
            out_file
                .write_all(format!("{}---{}\n", p5, p7).as_bytes())
                .unwrap()
        }
        // tend to the progress bar
        pb.inc(1);
        pb.set_message(&format!("{:3}%", 100 * i / length));
    }
    pb.finish_with_message(&"Done".to_string());
}

fn assemble_mapping(p5_path: &str, p7_path: &str, tsv_path: &str, dat_path: &str, out_path: &str) {
    let mapping = get_mapping(p5_path, p7_path);
    let clusters = StacksCluster::from_csv_records(tsv_path, dat_path);
    assemble_clusters(mapping, clusters, out_path);
}

fn count_locus_sizes(tsv_path: &str, dat_path: &str) {
    let _clusters = StacksCluster::from_csv_records(tsv_path, dat_path);
}

#[derive(StructOpt, Debug)]
#[structopt(name = "analyse stacks")]
enum StacksAnalyzer {
    GetMapping {
        /// Input files
        #[structopt(name = "p5-fq.gz-file")]
        p5_path: String,

        #[structopt(name = "p7-fq.gz-file")]
        p7_path: String,

        #[structopt(name = "tsv-file")]
        tsv_path: String,

        #[structopt(short = "d", long = "dat", required = true)]
        dat_path: String,

        #[structopt(short = "o", long = "out", required = true)]
        out_path: String,
    },

    CountLocusSizes {
        /// Input files
        #[structopt(name = "tsv-file")]
        tsv_path: String,

        #[structopt(short = "d", long = "dat", required = true)]
        dat_path: String,
    },
}

fn main() {
    match StacksAnalyzer::from_args() {
        StacksAnalyzer::GetMapping {
            p5_path,
            p7_path,
            tsv_path,
            dat_path,
            out_path,
        } => {
            assemble_mapping(&p5_path, &p7_path, &tsv_path, &dat_path, &out_path);
        }
        StacksAnalyzer::CountLocusSizes { tsv_path, dat_path } => {
            count_locus_sizes(&tsv_path, &dat_path);
        }
    }
}
