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

    #[arg(long = "report_reads")]
    get_read_num: bool,
}

fn main() {
    let args = Args::parse();
    let input = args.file;
    let sample_size = args.sample_size;

    let mut umis: IndexMap<i32, Vec<String>> = IndexMap::new();

    let outfile = format!("{}{}", input.rsplit_once('.').unwrap().0, "_umis.csv").to_string();
    let outfile = Path::new(&outfile);

    let _ = File::create(outfile);

    let mut num_reads: i64 = 0;

    if input.ends_with(".bam") {
        let bam_file = BamReader::from_path(&input, 8).unwrap();

        match args.sum {
            false => pull_umis_bam(bam_file, &mut umis, &args.separator, &mut num_reads),
            true => pull_umis_unsorted_bam(bam_file, &mut umis, &args.separator, &mut num_reads),
        }
    } else if input.ends_with(".txt") {
        pull_umis_txt(&Path::new(&input), &mut umis, &mut num_reads);
    }

    println!("Processing {} reads...", num_reads);
    match args.sample_size {
        0 => println!("Subsampling: None"),
        _ => println!(
            "Subsampling: {} read barcodes per position",
            args.sample_size
        ),
    };

    let process = match args.output {
        Output::Barcode => extract_umis,
        Output::Mean => extract_means,
        Output::Dist => extract_dist,
    };

    let position_reports = process(umis, sample_size);

    write_report(
        position_reports,
        outfile,
        args.get_read_num,
        args.sum,
        num_reads,
    );
}
