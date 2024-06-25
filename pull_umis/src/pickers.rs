use polars::frame::DataFrame;
use polars::lazy::dsl::mean;
use polars::prelude::*;
use rand::rngs::StdRng;
use rayon::prelude::*;
use std::fs::File;
use std::path::Path;
use std::sync::{Arc, Mutex};
use strsim::hamming;

use crate::{get_mean, utils::*};

use rand::seq::SliceRandom;
use rand::SeedableRng;

use indexmap::IndexMap;

// write the UMI barcodes themselves to file, with one column for each reference coordinate
pub fn extract_umis(mut store: IndexMap<i32, Vec<String>>, outfile: &Path, sample_size: usize) {
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

    write_report(dfs, outfile);
}

// write the mean edit distance to file. Mean is calculated per reference coordinate
// Note that 1000 UMIs will generate ~500_000 edit distance values, so using --sample is
// recommended
pub fn extract_means(mut store: IndexMap<i32, Vec<String>>, outfile: &Path, sample_size: usize) {
    let mut file = File::create(outfile).expect("Could not create file!");
    let mut dfs: Arc<Mutex<Vec<DataFrame>>> = Arc::new(Mutex::new(Vec::new()));

    store.par_drain(0..).for_each(|(pos, umi_list)| {
        let mut rng = StdRng::seed_from_u64(5);
        let mut edits: Vec<usize> = Vec::new();
        let sample = umi_list
            .choose_multiple(&mut rng, sample_size)
            .collect::<Vec<&String>>();

        let mut i = 0;
        for umi in &sample {
            for neighbor in &sample[i + 1..] {
                edits.push(hamming(&umi, neighbor).unwrap())
            }
            i += 1;
        }
        if !edits.is_empty() {
            let df = DataFrame::new(vec![Series::new(
                &pos.to_string(),
                [get_mean(edits) as f32],
            )])
            .unwrap();

            dfs.lock().unwrap().push(df);
        }
    });

    let dfs = Arc::try_unwrap(dfs).unwrap().into_inner().unwrap();

    write_report(dfs, outfile);
}

// write the frequency of edit distances 0-13 per reference coordinate.
// Note that 1000 UMIs will generate ~500_000 edit distance values, so using --sample is
// recommended
pub fn extract_dist(mut store: IndexMap<i32, Vec<String>>, outfile: &Path, sample_size: usize) {
    let mut file = File::create(outfile).expect("Could not create file!");
    let dfs: Arc<Mutex<Vec<DataFrame>>> = Arc::new(Mutex::new(Vec::new()));

    store.par_drain(0..).for_each(|(pos, umi_list)| {
        let mut rng = StdRng::seed_from_u64(5);
        let mut edits: Vec<usize> = Vec::new();
        let sample = umi_list
            .choose_multiple(&mut rng, sample_size)
            .collect::<Vec<&String>>();

        let mut i = 0;
        for umi in &sample {
            for neighbor in &sample[i + 1..] {
                edits.push(hamming(&umi, neighbor).unwrap())
            }
            i += 1;
        }

        let mut df = DataFrame::default();
        let _ = df.with_column(Series::new("edit", (0..13).collect::<Vec<i32>>()));

        let mut edit_col = (0..13).map(|_x| None).collect::<Vec<Option<i32>>>();

        if !edits.is_empty() {
            get_dist(&mut edits).drain(0..).for_each(|(edit, freq)| {
                edit_col[edit] = Some(freq);
            });

            let _ = df.with_column(Series::new(&pos.to_string(), edit_col));

            dfs.lock().unwrap().push(df);
        };
    });

    let mut first = true;
    dfs.lock().unwrap().iter_mut().for_each(|df| {
        if first {
            first = false;
        } else {
            let _ = df.drop_in_place("edit");
        }
    });

    let dfs = Arc::try_unwrap(dfs).unwrap().into_inner().unwrap();

    write_report(dfs, outfile);
}
