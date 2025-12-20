import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter raw counts")
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Input file path",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output file path",
    )
    parser.add_argument(
        "--sep",
        type=str,
        default="\t",
        help="Separator used in the input file",
    )

    args = parser.parse_args()

    # Read the input file
    df = pd.read_csv(args.input_file, sep=args.sep, index_col=0)

    # Remove genes that have count 0 in 95% or more of the samples
    df_filtered = df.loc[(df != 0).mean(axis=1) > 0.05]

    # Save the filtered DataFrame to a new file
    df_filtered.to_csv(args.output_file, sep=",")

    print(f"Original gene count: {df.shape[0]}")
    print(f"Filtered gene count: {df_filtered.shape[0]}")

    print(f"Original sample count: {df.shape[1]}")
    print(f"Filtered sample count (should be the same): {df_filtered.shape[1]}")

    # Check if there are samples with value 0 for all genes
    zero_samples = df_filtered.columns[(df_filtered.sum(axis=0) == 0)]

    if len(zero_samples) > 0:
        print("WARNING: The following samples have zero counts for all genes:")
        for s in zero_samples:
            print(f"  - {s}")
    else:
        print("OK: No samples with zero counts across all genes.")