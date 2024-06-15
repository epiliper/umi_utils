use rand::seq::SliceRandom;
use rand::{thread_rng, Rng};
use std::fs::{File, OpenOptions};
use std::mem;
use std::path::Path;
use triple_accel::hamming;
use indexmap::IndexMap;
use std::process;

use std::io::Write;

const BASES: [&str; 4] = [&"A", &"G", &"C", &"T"];

fn main() {
    let umis = generate_umis(8, 100);
    let amplicons = amplify_umis(umis, 8, 8, 1e-5, 0.8);
    println!("Number of amplicons: {}", amplicons.len());
    sequence_umis(amplicons, 8, 100, 1e-3);
}

pub fn generate_umis(umi_length: usize, num_umis: usize) -> Vec<String> {
    let mut rng = thread_rng();
    let mut umis: Vec<String> = Vec::with_capacity(num_umis);

    for num in 0..num_umis {
        let mut barcode = String::new();
        for bp in 0..umi_length {
            barcode.push_str(*BASES.choose(&mut rng).unwrap())
        }
        umis.push(barcode);
    }
    return umis;
}

pub fn amplify_umis(
    mut umis: Vec<String>,
    umi_length: usize,
    cycles: usize,
    error_rate: f32,
    amp_rate: f32,
) -> Vec<String> {
    let mut rng = thread_rng();

    for cycle in 0..cycles {
        let mut products: Vec<String> = Vec::new();

        for umi in &umis {
            if (rng.gen_range(0.0..1.0) as f32).to_bits() < amp_rate.to_bits() {
                let mut amp = umi.clone();

                for i in 0..umi_length {
                    if (rng.gen_range(0.0..1.0) as f32).to_bits() < error_rate.to_bits() {
                        print!("ERROR");
                        amp.replace_range(i..i + 1, *BASES.choose(&mut rng).unwrap());
                    }
                }
                products.push(amp);
            }
        }
        umis.append(&mut products);
    }
    return umis;
}

pub fn sequence_umis(mut umis: Vec<String>, umi_length: usize, depth: usize, error_rate: f32) {
    let mut reads: Vec<String> = Vec::new();

    umis.drain(0..).for_each(|umi| {
        let mut seqs: Vec<String> = Vec::new();
        let mut rng = thread_rng();

        for round in 0..depth {
            let mut read = umi.clone();
            for i in 0..umi_length {
                if (rng.gen_range(0.0..1.0) as f32).to_bits() < error_rate.to_bits() {
                    print!("error");
                    read.replace_range(i..i + 1, *BASES.choose(&mut rng).unwrap());
                }
            }
            seqs.push(read);
        }
        seqs.push(umi);
        reads.push(get_read(umi_length, seqs));
    });

    get_edit_distances(reads);
}

pub fn get_read(umi_length: usize, seqs: Vec<String>) -> String {

    let mut read = String::new();

    for i in 0..umi_length {

        let mut tally: IndexMap<char, i32> = IndexMap::new();

        for seq in &seqs {
            let hit = seq.as_bytes()[i] as char;
            tally.entry(hit).or_insert(0);
            *tally.get_mut(&hit).unwrap() += 1;
        }

        tally.sort_by(|a, b, c, d| d.cmp(b));
        let winning_base = tally.iter().next().unwrap().0;
        read.push(*winning_base);
    }
    return read
}

pub fn get_edit_distances(umis: Vec<String>) {
    let report_file = Path::new("edit_distance.txt");

    if !report_file.exists() {
        let _ = File::create(&report_file);
    }

    let mut f = OpenOptions::new()
        .write(true)
        .append(true)
        .open(&report_file)
        .expect("unable to open file");

    let mut i = 0;
    let bound = umis.len() - 1;

    while i < bound {
        for umi in &umis[0..] {
            for next_umi in &umis[i + 1..] {
                f.write(format!("{}\t\n", hamming(umi.as_bytes(), next_umi.as_bytes())).as_bytes())
                    .unwrap();
            }
            i += 1;
        }
    }
    mem::drop(umis);
}
