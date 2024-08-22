import pysam
from collections import Counter
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("bam")
parser.add_argument("--csv_name", default=None, required=False)
args = parser.parse_args()

TAG = "UG"


def count_reads_per_group(bam: str):
    tag_counts = Counter()
    with pysam.AlignmentFile(bam, "rb") as bam_file:
        for read in bam_file:
            if read.has_tag(TAG):
                tag_counts[read.get_tag(TAG)] += 1

    return tag_counts


# can use this data to make a histogram
def write_counts(bam_name, tag_counts: Counter, csv_name: str):
    if csv_name is None:
        csv_name = f"{os.path.join(
            os.path.dirname(bam_name), os.path.basename(bam_name).split(".")[0]
        )}_counts.csv"
    df = pd.DataFrame()
    df["counts"] = tag_counts.values()
    df.to_csv(csv_name, index=False)


write_counts(args.bam, count_reads_per_group(args.bam), args.csv_name)
