use bam::BamReader;
use clap::{Parser, ValueEnum};
use std::fs::File;
use std::path::Path;
mod pickers;
use pickers::*;
mod utils;
use indexmap::IndexMap;
use utils::*;

#[derive(ValueEnum, Debug, Clone)]
enum Output {
    Barcode,
    Dist,
    Mean,
}

#[derive(Parser, Debug)]
#[command(term_width = 0)]
struct Args {
    file: String,
    #[arg(long = "sep")]
    separator: String,

    #[arg(long = "sample", default_value_t = 0)]
    sample_size: usize,

    #[arg(long = "stat")]
    output: Output,
}

fn main() {
    let args = Args::parse();
    let input = args.file;
    let sample_size = args.sample_size;

    let bam_file = BamReader::from_path(&input, 8).unwrap();

    let outfile = format!("{}{}", input.split('.').next().unwrap(), "_umis.csv").to_string();
    let outfile = Path::new(&outfile);
    let _ = File::create(outfile);

    let mut umis: IndexMap<i32, Vec<String>> = IndexMap::new();

    for r in bam_file {
        let read = &r.unwrap();
        pull_umi(read, &mut umis, &args.separator)
    }

    let process = match args.output {
        Output::Barcode => extract_umis,
        Output::Mean => extract_means,
        Output::Dist => extract_dist,
    };

    process(umis, &outfile, sample_size);
}
