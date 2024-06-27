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
    #[arg(long = "sep", default_value = ":")]
    separator: String,

    #[arg(long = "sample", default_value_t = 0)]
    sample_size: usize,

    #[arg(long = "stat")]
    output: Output,

    #[arg(long = "sum")]
    sum: bool,
}

fn main() {
    let args = Args::parse();
    let input = args.file;
    let sample_size = args.sample_size;

    let mut umis: IndexMap<i32, Vec<String>> = IndexMap::new();

    let outfile = format!("{}{}", input.split('.').next().unwrap(), "_umis.csv").to_string();
    let outfile = Path::new(&outfile);

    let _ = File::create(outfile);

    if input.ends_with(".bam") {
        let bam_file = BamReader::from_path(&input, 8).unwrap();

        for r in bam_file {
            let read = &r.unwrap();
            pull_umi(read, &mut umis, &args.separator)
        }
    } else if input.ends_with(".txt") {
        pull_umis_txt(&Path::new(&input), &mut umis);
    }

    let process = match args.output {
        Output::Barcode => extract_umis,
        Output::Mean => extract_means,
        Output::Dist => extract_dist,
    };

    let position_reports = process(umis, &outfile, sample_size);

    write_report(position_reports, outfile, args.sum);
}
