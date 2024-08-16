use bam::BamReader;
use bam::Record;
use indexmap::IndexMap;
use polars::frame::DataFrame;
use polars::functions::concat_df_horizontal;
use polars::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// count the number of comparisons with edit distance of 0, 1, 2, etc.

pub fn get_read_pos(read: &Record) -> Option<i32> {
    if read.flag().is_reverse_strand() {
        let mut start = read.calculate_end();

        if !read.cigar().is_empty() {
            // set end pos as start to group with forward-reads covering same region
            start += read.cigar().soft_clipping(false) as i32; // pad with right-side soft clip
            return Some(start);
        }
    } else {
        let mut start = read.start();

        if !read.cigar().is_empty() {
            start -= read.cigar().soft_clipping(true) as i32; // pad with left-side soft clip
            return Some(start);
        }
    };

    return None;
}

pub fn get_dist(edits: &mut Vec<usize>) -> IndexMap<usize, i32> {
    let mut dist: IndexMap<usize, i32> = IndexMap::new();

    edits.drain(0..).for_each(|edit| {
        dist.entry(edit as usize).or_insert(0);

        *dist.get_mut(&edit).unwrap() += 1;
    });
    return dist;
}

// retrieve UMI from read record
fn get_umi(record: &Record, separator: &String) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(separator).last().unwrap().to_string()
}

// store UMI by position of associated read
pub fn pull_umis_bam(
    bam: BamReader<File>,
    store: &mut IndexMap<i32, Vec<String>>,
    separator: &String,
    num_reads: &mut i64,
) {
    for r in bam {
        let read = &r.unwrap();
        // if read is reverse to reference, group it by its last aligned base to the reference
        if read.flag().is_mapped() {
            let pos = get_read_pos(read).expect("ERROR: mapped read does not have CIGAR string!");
            store.entry(pos).or_insert(Vec::new());
            store.get_mut(&pos).unwrap().push(get_umi(read, separator));
            *num_reads += 1;
        }
    }
}

pub fn pull_umis_unsorted_bam(
    bam: BamReader<File>,
    store: &mut IndexMap<i32, Vec<String>>,
    separator: &String,
    num_reads: &mut i64,
) {
    for r in bam {
        let read = &r.unwrap();
        if read.flag().is_mapped() {
            store.entry(0).or_insert(Vec::new());
            store.get_mut(&0).unwrap().push(get_umi(read, separator));
            *num_reads += 1;
        }
    }
}

pub fn pull_umis_txt(txt: &Path, store: &mut IndexMap<i32, Vec<String>>, num_reads: &mut i64) {
    store.entry(1).or_insert(Vec::new());

    let file = File::open(txt).expect("unable to open txt file!");
    BufReader::new(file).lines().for_each(|line| {
        store.get_mut(&1).unwrap().push(line.unwrap().to_string());
        *num_reads += 1;
    });
}

// use only a subset of UMIs for downstream analysis, useful when working with high read depth
// recommended --sample value for BAM files >1GB: 5000
pub fn subsample(mut df: DataFrame, sample_size: usize) -> DataFrame {
    let col = &df.get_columns()[0];
    let len = col.len();

    if sample_size == 0 {
        return df;
    }

    if len >= sample_size {
        df = df
            .sample_n_literal(sample_size, false, false, Some(0))
            .unwrap();
        return df;
    }
    return df;
}

pub fn write_report(
    dfs: Vec<DataFrame>,
    outfile: &Path,
    get_reads: bool,
    summarize: bool,
    num_reads: i64,
) {
    let mut file = File::create(outfile).expect("Could not create file!");
    let mut new_report = concat_df_horizontal(&dfs).unwrap();
    new_report.align_chunks();

    if summarize {
        new_report = sum_cols(new_report, outfile);
    }

    if get_reads {
        new_report
            .with_column(Series::new("num_reads", [num_reads]))
            .unwrap();
    }

    CsvWriter::new(&mut file)
        .n_threads(8)
        .finish(&mut new_report)
        .unwrap();
}

// sum all edit distances across all recorded positions for genome-wide statistics
pub fn sum_cols(mut df: DataFrame, outfile: &Path) -> DataFrame {
    let mut col_names = df.get_column_names();
    col_names.remove(0); // remove "edit" column from summing

    let sample_name = outfile.to_str().unwrap().split("_umis.csv").next().unwrap();

    let sum = df
        .select(col_names)
        .unwrap()
        .sum_horizontal(polars::frame::NullStrategy::Ignore)
        .unwrap()
        .unwrap()
        .rename(sample_name)
        .clone()
        .fill_null(FillNullStrategy::Zero)
        .unwrap();

    let edit_col = df.drop_in_place("edit").unwrap();

    let summary = DataFrame::new(vec![edit_col, sum]).unwrap();
    return summary;
}

pub fn get_mean(list: Vec<usize>) -> f32 {
    return list.iter().sum::<usize>() as f32 / list.len() as f32;
}
