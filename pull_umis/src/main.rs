use bam::{BamReader, Record};
use clap::Parser;
use polars::frame::DataFrame;
use polars::functions::concat_df_horizontal;
use polars::prelude::*;
use rand::rngs::StdRng;
use rayon::prelude::*;
use std::fs::File;
use std::path::Path;
use strsim::hamming;
use std::sync::{Mutex, Arc};

use rand::seq::SliceRandom;
use rand::SeedableRng;

use indexmap::IndexMap;

//  const sample_size: usize = 1000;

#[derive(Parser, Debug)]
#[command(term_width = 0)]
struct Args {
    file: String,
    #[arg(long = "sep")]
    separator: String,

    #[arg(long = "sample")]
    sample_size: usize,
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

    // write_umis(umis, &outfile, sample_size);
    write_edit_distance(umis, &outfile, sample_size);
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

pub fn subsample(mut df: DataFrame, sample_size: usize) -> DataFrame {
    let col = &df.get_columns()[0];
    let len = col.len();

    if len >= sample_size {
        df = df
            .sample_n_literal(sample_size, false, false, Some(0))
            .unwrap();
        return df;
    }
    return df;
}

pub fn write_umis(mut store: IndexMap<i32, Vec<String>>, outfile: &Path, sample_size: usize) {
    let mut file = File::create(outfile).expect("Could not create file!");
    let mut dfs: Vec<DataFrame> = Vec::new();

    for (pos, umis) in store.drain(0..) {
        dfs.push(
            subsample(
                DataFrame::new(vec![Series::new(&pos.to_string(), umis)])
                    .unwrap()
                    .drop_nulls::<String>(None)
                    .unwrap(),
                sample_size, // This step saves ~50% file size, makes columns
            ), // disitinct lengths
        );
    }

    let mut new_report = concat_df_horizontal(&dfs).unwrap();
    new_report.align_chunks();

    CsvWriter::new(&mut file)
        .n_threads(8)
        .finish(&mut new_report)
        .unwrap();
}

pub fn write_edit_distance(
    mut store: IndexMap<i32, Vec<String>>,
    outfile: &Path,
    sample_size: usize,
) {


    let mut file = File::create(outfile).expect("Could not create file!");
    let mut dfs: Arc<Mutex<Vec<DataFrame>>> = Arc::new(Mutex::new(Vec::new()));

    store.par_drain(0..).for_each(|(pos, umi_list)| {
        let mut rng = StdRng::seed_from_u64(5);
        let mut edits: Vec<usize> = Vec::new();
        let sample = umi_list.choose_multiple(&mut rng, sample_size).collect::<Vec<&String>>();

        let mut i = 0;
        for umi in &sample {
            for neighbor in &sample[i+1..] {
                edits.push(hamming(&umi, neighbor).unwrap())
            }
            i += 1;
        }
        if !edits.is_empty() {
        let df = DataFrame::new(
            vec![Series::new(&pos.to_string(), [mean(edits) as f32])]
        ).unwrap();

        dfs.lock().unwrap().push(df);
        }
    });

    let mut new_report = concat_df_horizontal(&dfs.lock().unwrap()).unwrap();
    new_report.align_chunks();

    CsvWriter::new(&mut file)
        .n_threads(8)
        .finish(&mut new_report)
        .unwrap();
}

fn mean(list: Vec<usize>) -> f32 {

    return list.iter().sum::<usize>() as f32 / list.len() as f32
}
