use bam::{BamReader, Record};
use std::fs::{File, OpenOptions};
use std::path::Path;

use clap::Parser;
use std::io::Write;

use indexmap::IndexMap;

fn get_umi(record: &Record, separator: &String) -> String {
    let umi = String::from_utf8(record.name().to_vec());
    umi.unwrap().split(separator).last().unwrap().to_string()
}

pub fn pull_read(read: &Record, store: &mut IndexMap<i32, Vec<String>>, separator: &String) {
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

pub fn write_umis(position: i32, umis: Vec<String>, outfile: &Path) {

    let mut f = OpenOptions::new()
        .write(true)
        .append(true)
        .open(&outfile)
        .expect("unable to open file");

    f.write(&format!("{}\t", position).into_bytes());

    umis.iter().for_each(
        |umi| {
            f.write(&format!("{}\t", umi).into_bytes());
        });

    f.write(b"\n");
}

#[derive(Parser, Debug)]
#[command(term_width = 0)]
struct Args {
    file: String,
    #[arg(short = 's')]
    separator: String,
}

fn main() {
    let args = Args::parse();
    let input = args.file;

    let mut bam_file = BamReader::from_path(&input, 8).unwrap();

    let outfile = format!("{}{}", input.split('.').next().unwrap(), "_umis.tsv").to_string();
    let mut outfile = Path::new(&outfile);
    File::create(outfile);

    let mut umis: IndexMap<i32, Vec<String>> = IndexMap::new();

    for r in bam_file {
        let read = &r.unwrap();
        pull_read(read, &mut umis, &args.separator)
    }

    umis.drain(0..).for_each(
        |position| {
            write_umis(position.0, position.1, outfile)
        }
    )

}
