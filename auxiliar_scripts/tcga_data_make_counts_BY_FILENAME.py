import pandas as pd
import argparse
import os

def process_tcga_data(input_dir, df_sample_sheet, output_dir):
    # List the directories in the input directory
    dir_list = os.listdir(input_dir)
    dir_list = [d for d in dir_list if os.path.isdir(os.path.join(input_dir, d))]

    # DataFrames to hold counts
    df_rna_tumor = pd.DataFrame()
    df_rna_normal = pd.DataFrame()
    df_mirna_tumor = pd.DataFrame()
    df_mirna_normal = pd.DataFrame()

    df_rna_pooled = pd.DataFrame()
    df_mirna_pooled = pd.DataFrame()

    for dir in dir_list:
        # The dir will be listed in the File ID column of the sample sheet.
        # First Check if the dir is miRNA or RNA-Seq, get the row that dir comes from
        df_row = df_sample_sheet[df_sample_sheet["File ID"] == dir]
        condition = df_row["Tissue Type"].values[0]
        sample_name = df_row["Sample ID"].values[0]
        filename = df_row["File Name"].values[0]
        # Check if the row is RNA or miRNA
        if df_row["Data Type"].values[0] == "Gene Expression Quantification":
            # If it is RNA, then the count matrix file ends in TSV
            # and we have to ignore the first 6 lines of the file
            file_path = os.path.join(input_dir, dir, filename)
            df_counts = pd.read_csv(file_path, sep="\t", skiprows=6,
                                    names=["GeneID", "Symbol", "biotype", "unstranded",
                                    "stranded_first", "stranded_second", "tpm_unstranded",
                                    "fpkm_unstranded", "fpkm_uq_unstranded"])
            # Keep only protein coding genes
            df_counts = df_counts[df_counts["biotype"] == "protein_coding"]
            #print("RNA")
            #print(df_counts.head())
            #print("---------------------------------")

            # Add the counts to the RNA pooled dataframe
            # The RNA pooled dataframe will have GeneID as index, and the sampleID as columns
            df_counts = df_counts[["GeneID", "unstranded"]]
            df_counts = df_counts.set_index("GeneID")
            df_rna_pooled[filename] = df_counts["unstranded"]
            
            if condition == "Tumor":
                df_rna_tumor[filename] = df_counts["unstranded"]
            else:
                df_rna_normal[filename] = df_counts["unstranded"]

        else:
            # If it is miRNA, then the count matrix file ends in txt
            file_path = os.path.join(input_dir, dir, filename)
            df_counts = pd.read_csv(file_path, sep="\t")
            #print("miRNA")
            #print(df_counts.head())
            #print("---------------------------------")

            # Same logic as above, but for miRNA
            df_counts = df_counts[["miRNA_ID", "read_count"]]
            df_counts = df_counts.set_index("miRNA_ID")
            df_mirna_pooled[filename] = df_counts["read_count"]

            if condition == "Tumor":
                df_mirna_tumor[filename] = df_counts["read_count"]
            else:
                df_mirna_normal[filename] = df_counts["read_count"]

    # After processing all directories, save the dataframes to TSV files
    df_rna_tumor.to_csv(os.path.join(output_dir, "RNA_tumor_counts.tsv"), sep="\t")
    df_rna_normal.to_csv(os.path.join(output_dir, "RNA_normal_counts.tsv"), sep="\t")
    df_mirna_tumor.to_csv(os.path.join(output_dir, "miRNA_tumor_counts.tsv"), sep="\t")
    df_mirna_normal.to_csv(os.path.join(output_dir, "miRNA_normal_counts.tsv"), sep="\t")
    df_rna_pooled.to_csv(os.path.join(output_dir, "RNA_merged_counts.tsv"), sep="\t")
    df_mirna_pooled.to_csv(os.path.join(output_dir, "miRNA_merged_counts.tsv"), sep="\t")


def main():
    # Receive command line arguments
    parser = argparse.ArgumentParser(
        description="Process TCGA data files to create a counts matrix." \
        "The input directory should contain a list of other directories," \
        "each containing TCGA data files." \
        "The output will be 6 count matrices. RNA-tumor, RNA-normal," \
        "miRNA-tumor, miRNA-normal, RNA-merged, miRNA-merged.")
    
    parser.add_argument("input_dir", help="Directory containing TCGA directories")
    parser.add_argument("input_sample_sheet", help="Path to the TCGA sample sheet TSV file")
    parser.add_argument("output_dir", help="Directory to save the output count matrices")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Read the sample sheet
    df_sample_sheet = pd.read_csv(args.input_sample_sheet, sep="\t")

    # Process the TCGA data files
    process_tcga_data(args.input_dir, df_sample_sheet, args.output_dir)

if __name__ == "__main__":
    main()