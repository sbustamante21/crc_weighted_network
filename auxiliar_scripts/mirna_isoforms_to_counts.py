import argparse
import pandas as pd
import os
from pathlib import Path

def process_one_file(df):
    df = df[df['miRNA_region'].str.contains("MIMAT", na=False)].copy() # ONLY MIMAT rows
    mimats = []
    for el in df["miRNA_region"]:
        mimat = el.split(",")[1]
        if not mimat in mimats:
            mimats.append(mimat)
    mimats = sorted(mimats)

    counts = {}
    for mimat in mimats:
        # get the rows in which that mimat appears
        mini_df = df[df['miRNA_region'].str.contains(mimat, na=False)].copy()
        counts[mimat] = sum(mini_df["read_count"])

    return counts

def main():
    # Receive command line arguments
    parser = argparse.ArgumentParser(
        description="Process TCGA mirna-isoform data files to create a counts matrix."\
        "The resulting matrix will have MIMATS as row names."
    )
    parser.add_argument("input_dir", help="Directory containing TCGA isoforms count matrices")
    parser.add_argument("output_file", help="filename you will use for the output")
    args = parser.parse_args()
    
    dir_list = list(Path(args.input_dir).rglob("*quantification.txt"))

    all_counts = []
    for file in dir_list:
        df = pd.read_csv(file, sep="\t")
        counts = process_one_file(df)
        all_counts.append(counts)
    final_df = pd.DataFrame(all_counts).T
    final_df.columns = [file.name for file in dir_list]
    final_df = final_df.fillna(0)
    final_df = final_df.astype(int)

    final_df.to_csv(args.output_file, sep="\t", index=True, index_label="miRNA_ID")

if __name__ == "__main__":
    main()