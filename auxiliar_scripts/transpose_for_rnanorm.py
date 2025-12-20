# Transpose a TSV file, generate a CSV to later use in RNAnorm
import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transpose a TSV file for RNAnorm")
    parser.add_argument("--input_tsv", type=str, required=True,
                        help="Input TSV file to transpose")
    parser.add_argument("--sep", type=str, default="\t",
                        help="Separator used in the input TSV file")
    parser.add_argument("--output_csv", type=str, required=True,
                        help="Output CSV file after transposition")
    args = parser.parse_args()

    # Read the input TSV file
    df = pd.read_csv(args.input_tsv, sep=args.sep, index_col=0)

    # Transpose the DataFrame
    df_transposed = df.transpose()

    # Save the transposed DataFrame to CSV
    df_transposed.to_csv(args.output_csv)