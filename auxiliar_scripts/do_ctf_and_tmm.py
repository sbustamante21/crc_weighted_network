import pandas as pd
import argparse
from rnanorm import CTF
from rnanorm import TMM

# Receive a count matrix and return the CTF and TMM normalized versions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate CTF and TMM normalized count matrices," \
    "genes in rows, samples in columns")
    parser.add_argument("--counts", required=True, help="Path to the count matrix file.")
    parser.add_argument("--sep", default="\t", help="Separator used in the count matrix file.")
    parser.add_argument("--output_prefix", required=True, help="Prefix for the output files.")
    args = parser.parse_args()

    # Load count matrix
    counts_df = pd.read_csv(args.counts, sep=args.sep, index_col=0)
    counts_df = counts_df.T

    # add a pseudocount to avoid 0 values in the matrix
    counts_df = counts_df + 1

    # Calculate CTF and TMM normalized count matrices
    ctf_model = CTF(m_trim=0.3, a_trim=0.05).set_output(transform="pandas")
    tmm_model = TMM(m_trim=0.3, a_trim=0.05).set_output(transform="pandas")
    
    ctf_normalized = ctf_model.fit_transform(counts_df)
    tmm_normalized = tmm_model.fit_transform(counts_df)

    # Save the normalized count matrices
    ctf_normalized.to_csv(f"{args.output_prefix}_T_CTF.csv", sep=args.sep)
    tmm_normalized.to_csv(f"{args.output_prefix}_T_TMM.csv", sep=args.sep)
