import pandas as pd

# Load your metadata
metadata = pd.read_csv("your_metadata.csv", sep="\t")  # Adjust separator if needed

# -------------------------------
# 🔍 Filter for mRNA or miRNA
# -------------------------------
# Example: filter for mRNA
rna_meta = metadata[metadata['Data Type'].str.contains("Gene Expression", case=False)]

# Example: filter for miRNA
mirna_meta = metadata[metadata['Data Type'].str.contains("miRNA", case=False)]

# -------------------------------
# ✅ WGCNA: all samples (tumor or normal)
# -------------------------------
wgcna_rna_samples = rna_meta['Sample ID'].unique().tolist()
wgcna_mirna_samples = mirna_meta['Sample ID'].unique().tolist()

# -------------------------------
# ✅ DE: paired samples (same Case ID with both tumor and normal)
# -------------------------------
def get_paired_samples(df):
    # Count tissue types per Case ID
    tissue_counts = df.groupby(['Case ID', 'Tissue Type']).size().unstack(fill_value=0)
    paired_cases = tissue_counts[(tissue_counts.get('Tumor', 0) > 0) & (tissue_counts.get('Normal', 0) > 0)].index
    # Filter samples from paired cases
    paired_df = df[df['Case ID'].isin(paired_cases)]
    return paired_df['Sample ID'].unique().tolist()

de_rna_samples = get_paired_samples(rna_meta)
de_mirna_samples = get_paired_samples(mirna_meta)

# -------------------------------
# 🧾 Output summary
# -------------------------------
print("✅ WGCNA - RNA samples:", len(wgcna_rna_samples))
print("✅ WGCNA - miRNA samples:", len(wgcna_mirna_samples))
print("🧬 DE - Paired RNA samples:", len(de_rna_samples))
print("🧬 DE - Paired miRNA samples:", len(de_mirna_samples))