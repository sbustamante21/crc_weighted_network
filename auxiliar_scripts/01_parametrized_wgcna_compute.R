#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(WGCNA)
})

setwd("/home/sbustamante/WGCNA_colorectal_cancer/network_white_tcga")

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads=50)

# INPUTS
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    expr_file   <- args[2] # --expr
    traits_file <- args[4] # --traits
    outdir      <- args[6] # --outdir
}

# 1) Expression matrix file (samples x genes). Recommended: CSV with header + rownames as first column.
#expr_file   <- "count_matrices/sampledRNA_TMM_T_merged_counts.csv"
# 2) Traits file (samples x traits). Row names must match sample IDs in datExpr.
#traits_file <- "count_matrices/gdc_sample_sheet.2026-01-21_GENES+MIRNAISOFORMS.csv"   # can be NULL if you don't have traits yet
# 3) Output directory
#outdir <- "sampled_rna_pooled_tmm"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 4) Network settings
networkType <- "signed"
TOMType     <- "signed"
corType     <- "bicor"   # robust for RNA-seq-ish data
minModuleSize <- 30
mergeCutHeight <- 0.25

# 5) Soft-threshold search range
powers <- 1:20

# 6) Choose whether to save TOMs (big files)
saveTOMs <- TRUE

# 7) Cytoscape export thresholds
cytoscape_TOM_threshold <- 0.05   # lower => more edges, bigger files
top_genes_for_TOMplot   <- 500    # subset size for optional TOMplot data export

# =========================
# LOAD DATA
# =========================
message("Loading expression: ", expr_file)
datExpr <- read.table(expr_file, row.names = 1, sep=",", header=TRUE, check.names = FALSE)
# Ensure numeric matrix
#datExpr <- as.data.frame(datExpr)
#datExpr[] <- lapply(datExpr, function(x) as.numeric(as.character(x)))
#datExpr <- as.data.frame(datExpr)
#rownames(datExpr) <- rownames(read.table(expr_file, header=TRUE, sep=",", row.names=1, check.names=FALSE))

# Basic QC: good samples/genes
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  message("Removing bad samples/genes...")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Load traits if provided
datTraits <- NULL
if (!is.null(traits_file) && file.exists(traits_file)) {
  message("Loading traits: ", traits_file)
  #datTraits <- read.csv(traits_file, header=TRUE, sep=",", row.names = 1, check.names=FALSE) # TAIWAN
  datTraits <- read.csv(traits_file, header=TRUE, sep=",", row.names = 2, check.names=FALSE) # TCGA
  
  # Match samples
  common <- intersect(rownames(datExpr), rownames(datTraits))
  datExpr   <- datExpr[common, , drop=FALSE]
  datTraits <- datTraits[common, , drop=FALSE]
}

# Save cleaned inputs
saveRDS(datExpr, file = file.path(outdir, "datExpr_clean.rds"))
if (!is.null(datTraits)) saveRDS(datTraits, file = file.path(outdir, "datTraits_clean.rds"))

# =========================
# SOFT THRESHOLD
# =========================
message("Running pickSoftThreshold...")
sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         networkType = networkType,
                         corFnc = if (corType == "bicor") "bicor" else "cor",
                         verbose = 5)

saveRDS(sft, file = file.path(outdir, "softThreshold_sft.rds"))
write.table(sft$fitIndices,
            file = file.path(outdir, "softThreshold_fitIndices.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Choose a power automatically (simple heuristic): first power with R^2 >= 0.8, else best R^2
fit <- sft$fitIndices
signedR2 <- -sign(fit[,3]) * fit[,2]
candidate <- fit[which(signedR2 >= 0.8), 1]
softPower <- if (length(candidate) > 0) candidate[1] else fit[which.max(signedR2), 1]
message("Chosen softPower = ", softPower)
writeLines(paste0("softPower\t", softPower), con = file.path(outdir, "chosen_softPower.tsv"))

# =========================
# BLOCKWISE MODULES
# =========================
message("Running blockwiseModules (this is the heavy step)...")

net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = networkType,
  TOMType = TOMType,
  corType = corType,
  maxPOutliers = if (corType == "bicor") 0.1 else 1,
  minModuleSize = minModuleSize,
  reassignThreshold = 0,
  mergeCutHeight = mergeCutHeight,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = saveTOMs,
  saveTOMFileBase = file.path(outdir, "TOM"),
  verbose = 3
)

saveRDS(net, file = file.path(outdir, "blockwise_net.rds"))

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
MEs <- net$MEs

saveRDS(moduleColors, file = file.path(outdir, "moduleColors.rds"))
saveRDS(MEs,          file = file.path(outdir, "MEs.rds"))

# =========================
# MODULE GENE LISTS
# =========================
geneNames <- colnames(datExpr)

module_table <- data.frame(
  gene = geneNames,
  module = moduleColors,
  stringsAsFactors = FALSE
)

write.table(module_table,
            file = file.path(outdir, "module_membership_gene2module.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Per-module gene lists
dir.create(file.path(outdir, "modules"), showWarnings = FALSE)
for (m in sort(unique(moduleColors))) {
  genes_m <- module_table$gene[module_table$module == m]
  writeLines(genes_m, con = file.path(outdir, "modules", paste0("genes_", m, ".txt")))
}

nSamples <- nrow(datExpr)

# =========================
# kME (MODULE MEMBERSHIP) and P VALUES
# =========================
message("Computing kME (module membership)...")
kME <- as.data.frame(cor(datExpr, MEs, use="p"))
colnames(kME) <- paste0("kME_", colnames(MEs))

kME_P <- as.data.frame(corPvalueStudent(as.matrix(kME), nSamples=nrow(datExpr)))
colnames(kME_P) <- paste0(colnames(kME), "_P")

kME$gene <- geneNames
kME_P$gene <- geneNames

write.table(kME, file = file.path(outdir, "kME_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(kME, file = file.path(outdir, "kME_table.rds"))

write.table(kME_P, file= file.path(outdir, "kME_P_table.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(kME_P, file=file.path(outdir, "kME_P_table.rds"))

# =========================
# GENE SIGNIFICANCE and P values (optional, if traits exist)
# =========================
if (!is.null(datTraits)) {
  message("Computing Gene Significance (GS) vs traits...")
  # If traits are not numeric, convert factors to numeric as needed
  traits_num <- datTraits
  for (cn in colnames(traits_num)) {
    if (!is.numeric(traits_num[[cn]])) traits_num[[cn]] <- as.numeric(as.factor(traits_num[[cn]]))
  }

  GS <- as.data.frame(cor(datExpr, traits_num, use="p"))
  colnames(GS) <- paste0("GS_", colnames(traits_num))
  
  GS_P <- as.data.frame(corPvalueStudent(as.matrix(GS), nSamples=nSamples))
  colnames(GS_P) <- paste0(colnames(GS), "_P")

  GS$gene <- geneNames
  GS_P$gene <- geneNames

  write.table(GS, file = file.path(outdir, "GS_table.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(GS, file = file.path(outdir, "GS_table.rds"))

  write.table(GS_P, file = file.path(outdir, "GS_P_table.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(GS_P, file = file.path(outdir, "GS_P_table.rds"))

}

# =========================
# CYTOSCAPE EXPORT (per module)
# =========================
message("Preparing Cytoscape exports...")
dir.create(file.path(outdir, "cytoscape"), showWarnings = FALSE)

# Load TOM (can be huge). Only feasible if saveTOMs=TRUE
# We'll export for each module a subnetwork using TOM similarity.
if (saveTOMs) {
  # NOTE: blockwiseModules saves TOM per block. If you have multiple blocks,
  # full TOM reconstruction is non-trivial. For Cytoscape you can instead use
  # intramodular edges based on adjacency. We'll do adjacency-based export (lightweight).
  adj <- adjacency(datExpr, power = softPower, type = networkType, corFnc = if (corType == "bicor") "bicor" else "cor")
  for (m in setdiff(unique(moduleColors), "grey")) {
    inMod <- moduleColors == m
    modGenes <- geneNames[inMod]

    if (length(modGenes) < 5) next
    subAdj <- adj[inMod, inMod]
    colnames(subAdj) <- rownames(subAdj) <- modGenes

    cy <- exportNetworkToCytoscape(
      subAdj,
      edgeFile = file.path(outdir, "cytoscape", paste0("edges_", m, ".txt")),
      nodeFile = file.path(outdir, "cytoscape", paste0("nodes_", m, ".txt")),
      weighted = TRUE,
      threshold = cytoscape_TOM_threshold,
      nodeNames = modGenes,
      nodeAttr = moduleColors[inMod]
    )
  }
}

# =========================
# OPTIONAL: export a subset for TOMplot later (PC script)
# =========================
# Save a subset of genes for your PC to generate a TOM heatmap if you want.
set.seed(1)
subsetGenes <- geneNames[order(apply(datExpr, 2, var), decreasing = TRUE)][1:min(top_genes_for_TOMplot, ncol(datExpr))]
writeLines(subsetGenes, con = file.path(outdir, "TOMplot_subset_genes.txt"))

message("DONE. Outputs in: ", outdir)
