# This script takes the sample sheet and a directory and created the miRNA
# and RNA-Seq counts in different files

import pandas as pd
import argparse
import os

def main():
    # Receive the sample sheet file
    parser = argparse.ArgumentParser(description="Create RNA-Seq and miRNA-Seq counts")
    parser.add_argument("sample_sheet", help="Path to the sample sheet file")
    parser.add_argument("input_dir", help="Directory containing the count files, must be a directory full of directories")
    parser.add_argument("output_dir", help="Directory to save the output files")
    args = parser.parse_args()

    # create dataframe from sample sheet file
    samples_df = pd.read_csv(args.sample_sheet, sep="\t")

    # Create output directories if they don't exist
    os.makedirs(args.output_dir, exist_ok=True)


    # create empty RNA and miRNA count dataframes
    rna_counts = pd.DataFrame()
    mirna_counts = pd.DataFrame()

    for file_id in samples_df["File ID"]:
        # Get the sample ID and the file name for the current sample
        sample_id = samples_df.loc[samples_df["File ID"] == file_id, "Sample ID"].values[0]
        file_name = samples_df.loc[samples_df["File ID"] == file_id, "File Name"].values[0]

        # Build the path to the whole file, using file_id and file_name
        file_path = os.path.join(args.input_dir, file_id, file_name)

        # check if the row corresponding to that file_id has miRNA in the Data Type column
        if "miRNA" in samples_df.loc[samples_df["File ID"] == file_id, "Data Type"].values[0]:
            mirna_df = pd.read_csv(file_path, sep="\t")
            mirna_counts[sample_id] = mirna_df["read_count"]
        else:
            # if it is rna, ignore the first 6 lines
            rna_df = pd.read_csv(file_path, sep="\t", skiprows=6, names=
                                 ["gene_id", "gene_name", "gene_type", "unstranded",
                                  "stranded_first", "stranded_second",
                                  "tpm_unstranded", "fpkm_unstranded",
                                  "fpkm_uq_unstranded"])
            rna_df_protein_coding = rna_df[rna_df["gene_type"] == "protein_coding"]
            rna_counts[sample_id] = rna_df_protein_coding["unstranded"]
            
    mirna_counts.set_index(mirna_df["miRNA_ID"], inplace=True)
    rna_counts.set_index(rna_df_protein_coding["gene_name"], inplace=True)

    mirna_counts.to_csv(os.path.join(args.output_dir, "mirna_counts.tsv"), sep="\t")
    rna_counts.to_csv(os.path.join(args.output_dir, "rna_counts.tsv"), sep="\t")

if __name__ == "__main__":
    main()