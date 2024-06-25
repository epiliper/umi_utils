use bam::{BamReader, Record};
use indexmap::IndexMap;
use polars::frame::DataFrame;
use polars::functions::concat_df_horizontal;
use polars::prelude::*;
use std::fs::File;
use std::path::Path;

pub fn get_dist(edits: &mut Vec<usize>) -> IndexMap<usize, i32> {
    let mut dist: IndexMap<usize, i32> = IndexMap::new();

    edits.drain(0..).for_each(|edit| {
        dist.entry(edit as usize).or_insert(0);

        *dist.get_mut(&edit).unwrap() += 1;
    });
    return dist;
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

pub fn write_report(dfs: Vec<DataFrame>, outfile: &Path) {
    let mut file = File::create(outfile).expect("Could not create file!");
    let mut new_report = concat_df_horizontal(&dfs).unwrap();
    new_report.align_chunks();

    CsvWriter::new(&mut file)
        .n_threads(8)
        .finish(&mut new_report)
        .unwrap();
}

pub fn get_mean(list: Vec<usize>) -> f32 {
    return list.iter().sum::<usize>() as f32 / list.len() as f32;
}
