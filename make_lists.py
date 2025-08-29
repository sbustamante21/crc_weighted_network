import pandas as pd

# load metadata (tsv)
meta = pd.read_csv("metadata.tsv", sep="\t", dtype=str)

# normalize strings (optional)
meta["Data Type"] = meta["Data Type"].str.strip()
meta["Tissue Type"] = meta["Tissue Type"].str.strip()
meta["Case ID"] = meta["Case ID"].str.strip()
meta["Sample ID"] = meta["Sample ID"].str.strip()

# Map assay and condition
assay_map = {
    "Gene Expression Quantification": "RNA",
    "miRNA Expression Quantification": "miRNA"
}
cond_map = {
    "Primary Tumor": "Tumor",
    "Solid Tissue Normal": "Normal"
}

meta = meta[meta["Data Type"].isin(assay_map.keys())].copy()
meta["Assay"] = meta["Data Type"].map(assay_map)
meta["Condition"] = meta["Tissue Type"].map(cond_map)

# Remove rows with missing mapping
meta = meta[meta["Assay"].notna() & meta["Condition"].notna()]

# For convenience, create columns
meta["Case_Assay_Cond"] = meta["Case ID"] + "|" + meta["Assay"] + "|" + meta["Condition"]

# 1) Patients with RNA paired (both tumor & normal RNA)
rna = meta[meta["Assay"] == "RNA"]
rna_counts = rna.groupby("Case ID")["Condition"].unique().apply(list)
rna_paired_cases = [case for case, conds in rna_counts.items() if set(["Tumor", "Normal"]).issubset(set(conds))]
print("RNA paired cases (have both Tumor & Normal RNA):", len(rna_paired_cases))

# RNA paired sample IDs (both tumor and normal)
rna_paired_samples = meta[(meta["Case ID"].isin(rna_paired_cases)) & (meta["Assay"]=="RNA")]
# Optional: pivot to get sample IDs per case/condition
rna_pivot = rna_paired_samples.pivot_table(index="Case ID", columns="Condition", values="Sample ID", aggfunc=lambda x: ";".join(x))
rna_pivot.to_csv("rna_paired_sampleIDs_per_case.tsv", sep="\t")

# 2) Patients with miRNA paired
mirna = meta[meta["Assay"] == "miRNA"]
mirna_counts = mirna.groupby("Case ID")["Condition"].unique().apply(list)
mirna_paired_cases = [case for case, conds in mirna_counts.items() if set(["Tumor", "Normal"]).issubset(set(conds))]
print("miRNA paired cases:", len(mirna_paired_cases))

mirna_paired_samples = meta[(meta["Case ID"].isin(mirna_paired_cases)) & (meta["Assay"]=="miRNA")]
mirna_pivot = mirna_paired_samples.pivot_table(index="Case ID", columns="Condition", values="Sample ID", aggfunc=lambda x: ";".join(x))
mirna_pivot.to_csv("mirna_paired_sampleIDs_per_case.tsv", sep="\t")

# 3) Fully paired across both assays (RNA+miRNA, both Tumor & Normal)
# require presence of (RNA, Tumor), (RNA, Normal), (miRNA, Tumor), (miRNA, Normal)
cases = meta["Case ID"].unique().tolist()
fully_paired = []
for c in cases:
    sub = meta[meta["Case ID"] == c]
    has_rna_tumor = ((sub["Assay"]=="RNA") & (sub["Condition"]=="Tumor")).any()
    has_rna_normal = ((sub["Assay"]=="RNA") & (sub["Condition"]=="Normal")).any()
    has_mi_tumor = ((sub["Assay"]=="miRNA") & (sub["Condition"]=="Tumor")).any()
    has_mi_normal = ((sub["Assay"]=="miRNA") & (sub["Condition"]=="Normal")).any()
    if has_rna_tumor and has_rna_normal and has_mi_tumor and has_mi_normal:
        fully_paired.append(c)

print("Fully paired cases (RNA+miRNA, both Tumor & Normal):", len(fully_paired))

# 4) Samples usable for WGCNA (all RNA samples; all miRNA samples)
rna_wgcna_samples = meta[meta["Assay"] == "RNA"]["Sample ID"].unique().tolist()
mirna_wgcna_samples = meta[meta["Assay"] == "miRNA"]["Sample ID"].unique().tolist()
print("RNA WGCNA samples:", len(rna_wgcna_samples))
print("miRNA WGCNA samples:", len(mirna_wgcna_samples))

# 5) Samples usable for DE (paired): sample IDs per case to build paired DE datasets
# mRNA DE: keep sample IDs for cases in rna_paired_cases
rna_de_samples = meta[(meta["Case ID"].isin(rna_paired_cases)) & (meta["Assay"]=="RNA")]["Sample ID"].tolist()
# miRNA DE: keep sample IDs for cases in mirna_paired_cases
mirna_de_samples = meta[(meta["Case ID"].isin(mirna_paired_cases)) & (meta["Assay"]=="miRNA")]["Sample ID"].tolist()

# 6) Save results
pd.Series(rna_paired_cases).to_csv("rna_paired_cases.txt", index=False, header=False)
pd.Series(mirna_paired_cases).to_csv("mirna_paired_cases.txt", index=False, header=False)
pd.Series(fully_paired).to_csv("fully_paired_cases.txt", index=False, header=False)

pd.Series(rna_de_samples).to_csv("rna_de_sampleIDs.txt", index=False, header=False)
pd.Series(mirna_de_samples).to_csv("mirna_de_sampleIDs.txt", index=False, header=False)

pd.Series(rna_wgcna_samples).to_csv("rna_wgcna_sampleIDs.txt", index=False, header=False)
pd.Series(mirna_wgcna_samples).to_csv("mirna_wgcna_sampleIDs.txt", index=False, header=False)

print("Done. Files written: rna_de_sampleIDs.txt, mirna_de_sampleIDs.txt, rna_wgcna_sampleIDs.txt, mirna_wgcna_sampleIDs.txt, and paired case lists.")