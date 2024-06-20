use bam::{BamReader, Record};
use clap::Parser;
use polars::functions::concat_df_horizontal;
use polars::prelude::*;
use std::fs::File;
use std::path::Path;

use indexmap::IndexMap;

fn main() {
    let args = Args::parse();
    let input = args.file;

    let bam_file = BamReader::from_path(&input, 8).unwrap();

    let outfile = format!("{}{}", input.split('.').next().unwrap(), "_umis.csv").to_string();
    let outfile = Path::new(&outfile);
    let _ = File::create(outfile);

    let mut umis: IndexMap<i32, Vec<String>> = IndexMap::new();

    for r in bam_file {
        let read = &r.unwrap();
        pull_umi(read, &mut umis, &args.separator)
    }

    write_tsv(umis, &outfile);
}

fn get_umi(record: &Record, separator: &String) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(separator).last().unwrap().to_string()
}

pub fn pull_umi(read: &Record, store: &mut IndexMap<i32, Vec<String>>, separator: &String) {
    // if read is reverse to reference, group it by its last aligned base to the reference
    if read.flag().is_mapped() && read.flag().is_reverse_strand() {
        let pos = read.calculate_end() + 1;
        store.entry(pos).or_insert(Vec::new());
        store.get_mut(&pos).unwrap().push(get_umi(read, separator))
    }
    // otherwise, use its first position to reference
    else if read.flag().is_mapped() {
        let pos = read.start() + 1;
        store.entry(pos).or_insert(Vec::new());
        store.get_mut(&pos).unwrap().push(get_umi(read, separator))
    }
}

pub fn write_tsv(store: IndexMap<i32, Vec<String>>, outfile: &Path) {
    let mut file = File::create(outfile).expect("Could not create file!");
    let mut dfs: Vec<DataFrame> = Vec::new();

    for (pos, umis) in store {
        dfs.push(DataFrame::new(vec![Series::new(&pos.to_string(), umis)]).unwrap());
    }

    let mut new_report = concat_df_horizontal(&dfs).unwrap();
    new_report.align_chunks();

    CsvWriter::new(&mut file)
        .n_threads(8)
        .finish(&mut new_report)
        .unwrap();
}

#[derive(Parser, Debug)]
#[command(term_width = 0)]
struct Args {
    file: String,
    #[arg(short = 's')]
    separator: String,
}
