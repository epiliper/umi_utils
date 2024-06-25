use bam::{BamReader, Record};
use clap::{Parser, ValueEnum};
use polars::frame::DataFrame;
use polars::functions::concat_df_horizontal;
use polars::prelude::*;
use rand::rngs::StdRng;
use rayon::prelude::*;
use std::fs::File;
use std::path::Path;
use std::sync::{Arc, Mutex};
use strsim::hamming;

mod pickers;
use pickers::*;

mod utils;
use utils::*;

use rand::seq::SliceRandom;
use rand::SeedableRng;

use indexmap::IndexMap;

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

    #[arg(long = "sample")]
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

    // write_umis(umis, &outfile, sample_size);
    process(umis, &outfile, sample_size);
}
