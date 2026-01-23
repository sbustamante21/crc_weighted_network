import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Divide count matrices into tumor and normal based on metadata.")
    parser.add_argument("--counts", required=True, help="Path to the count matrix file. GENES ON COLUMNS")
    parser.add_argument("--sep", default=",", help="Separator used in the count matrix file.")
    parser.add_argument("--metadata", required=True, help="Path to the metadata file (CSV format).")
    parser.add_argument("--output_prefix", required=True, help="Prefix for the output files.")
    args = parser.parse_args()

    # Load count matrix
    counts_df = pd.read_csv(args.counts, sep=args.sep, index_col=0)
    counts_df = counts_df.T # because we receive genes on columns, here we change to samples on columns

    # Load metadata
    #metadata_df = pd.read_csv(args.metadata, index_col=0) # for TAIWAN data
    metadata_df = pd.read_csv(args.metadata, index_col=1) # for TCGA data

    # Ensure the order of samples in metadata matches the count matrix
    metadata_df = metadata_df.loc[counts_df.columns]

    # Divide samples into tumor and normal
    #tumor_samples = metadata_df[metadata_df['PHENOTYPE'] == 'neoplastic'].index # TAIWAN
    #normal_samples = metadata_df[metadata_df['PHENOTYPE'] == 'adjacent normal'].index # TAIWAN

    tumor_samples = metadata_df[metadata_df['Tissue.Type'] == 'Tumor'].index # TCGA
    normal_samples = metadata_df[metadata_df['Tissue.Type'] == 'Normal'].index # TCGA

    # Create tumor and normal count matrices
    tumor_counts_df = counts_df[tumor_samples]
    normal_counts_df = counts_df[normal_samples]

    # Save the divided count matrices, also the same input but transposed too
    counts_df.T.to_csv(f"{args.output_prefix}_T_merged_counts.csv", sep=args.sep)
    tumor_counts_df.T.to_csv(f"{args.output_prefix}_T_tumor_counts.csv", sep=args.sep)
    normal_counts_df.T.to_csv(f"{args.output_prefix}_T_normal_counts.csv", sep=args.sep)
