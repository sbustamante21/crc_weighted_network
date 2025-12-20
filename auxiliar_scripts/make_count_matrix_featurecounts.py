import argparse
import os
from os.path import isfile, join
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Generate count matrix from featureCounts output")
    parser.add_argument(
        "root_dir",
        help="Root directory containing featureCounts output for each sample",
    )
    parser.add_argument(
        "out_file",
        help="Name of the output file to generate",
    )
    args = parser.parse_args()

    files = [f for f in os.listdir(args.root_dir) if isfile(join(args.root_dir, f))]
    # read first file as tsv, ignore first line
    first_file_df = pd.read_csv(join(args.root_dir, files[0]), sep="\t", skiprows=[0])

    # take just the first column, name column "GeneID"
    first_file_df = first_file_df.iloc[:, 0]
    first_file_df.name = "GeneID"

    for f in files:
        # take last column from each file, column name will be first 11 characters from file name
        #sample_name = f[0: 10]
        sample_name = f.split("_")[0]

        df = pd.read_csv(join(args.root_dir, f), sep="\t", skiprows=[0])
        df = df.iloc[:, 6]
        df.name = sample_name
        first_file_df = pd.concat([first_file_df, df], axis=1)

    first_file_df.to_csv(args.out_file, index=False, sep="\t")

if __name__ == "__main__":
    main()