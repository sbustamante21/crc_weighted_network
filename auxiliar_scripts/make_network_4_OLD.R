library("DESeq2")
library("WGCNA")
library("pheatmap")
library("dynamicTreeCut")
library("grid")

# Load the counts and the metadata
counts <- read.table("../input_data/rna_counts.txt", header = TRUE, row.names = 1, sep="\t")
metadata <- read.table("../input_data/all_metadata.csv", header = TRUE, row.names = 1, sep=",")

combined_counts <- counts
zero_proportion <- rowMeans(combined_counts == 0)

# Keep only those rows where fewer than 95% of values are zero
combined_counts <- combined_counts[zero_proportion < 0.95, ]

# Metadata rownames must be combined counts colnames, and must be in the same order
metadata <- metadata[match(colnames(combined_counts), rownames(metadata)), ]

print("Performing differential expression analysis...")
dds <- DESeqDataSetFromMatrix(countData = combined_counts, # Differential Expression is made on the combined counts
                              colData = metadata,
                              design = ~ PHENOTYPE)
dds <- DESeq(dds)

res = results(dds, contrast = c("PHENOTYPE", "adjacent normal", "neoplastic")) # Modifiable based on the metadata
res = res[order(res$padj), ]
res <- res[!is.na(res$padj), ]

res_ordered <- res[order(res$padj), ]
df = res_ordered
df$Category <- "Not significant" # Criteria can be changed
df[which(df$log2FoldChange >= 1 & df$padj <= 0.001),]$Category <- "Up-regulated"
df[which(df$log2FoldChange <= -1 & df$padj <= 0.001),]$Category <- "Down-regulated"

# save the differential expression dataframe
write.csv(as.data.frame(df), "../output/differential_expression_results.csv")

df_sig <- df[df$Category != "Not significant", ]

sig_genes <- rownames(df_sig)

#norm_counts <- counts(dds, normalized = TRUE)[sig_genes, ] # THIS TO USE ONLY DEGs
norm_counts <- counts(dds, normalized=TRUE) # USE THIS TO CREATE A FULL NETWORK

# INICIO PARTE TEST
# Identify DEGs (e.g., up/down-regulated)
deg_genes <- rownames(df[df$Category != "Not significant", ])
# Identify non-DEGs
non_deg_genes <- setdiff(rownames(df), deg_genes)
# Randomly select a subset
set.seed(123)  # for reproducibility
deg_subset <- sample(deg_genes, size = 500)
non_deg_subset <- sample(non_deg_genes, size = 500)
# Combine
selected_genes <- c(deg_subset, non_deg_subset)
# Subset norm_counts
subset_counts <- norm_counts[selected_genes, ]
# FIN PARTE TEST

#datExpr0 <- t(subset_counts) # IF ONLY USING DEGs

datExpr0 <- t(norm_counts)

sampleTree <- hclust(dist(datExpr0), method = "average")

print("Sample clustering...")
svg("../output/sample_clustering.svg", width=8, height=6)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#Plot a line showing the cut-off
abline(h = 4.5e+06, col = "red") # CutHeight can vary
dev.off()

clust = cutreeStatic(sampleTree, cutHeight = 4.5e+06, minSize = 10)

# Cluster 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

stopifnot(all(rownames(datExpr0) %in% rownames(metadata)))

# Start from the columns we need (adjust if names differ)
datTraits <- metadata[, c("PHENOTYPE", "age", "disease_stage", "sex", "tissue"), drop = FALSE]

# Helpers
normlower <- function(x) tolower(trimws(as.character(x)))

# PHENOTYPE: adjacent normal -> 0, neoplastic -> 1
if ("PHENOTYPE" %in% names(datTraits)) {
  ph <- normlower(datTraits$PHENOTYPE)
  datTraits$PHENOTYPE <- ifelse(ph %in% c("adjacent normal","adjacent_normal","normal","adjacent"), 0,
                          ifelse(ph %in% c("neoplastic","tumor","tumour","cancer","tumoral"), 1, NA_real_))
}

# age: numeric
if ("age" %in% names(datTraits)) {
  datTraits$age <- suppressWarnings(as.numeric(datTraits$age))
}

# disease_stage: TNM1..TNM4 -> 1..4
if ("disease_stage" %in% names(datTraits)) {
  ds <- normlower(datTraits$disease_stage)
  # accept formats like "tnm1", "TNM2", " tnm3 "
  datTraits$disease_stage_num <- dplyr::case_when(
    grepl("tnm4|\\b4\\b", ds) ~ 4,
    grepl("tnm3|\\b3\\b", ds) ~ 3,
    grepl("tnm2|\\b2\\b", ds) ~ 2,
    grepl("tnm1|\\b1\\b", ds) ~ 1,
    TRUE ~ NA_real_
  )
}

# sex: male -> 1, female -> 0
if ("sex" %in% names(datTraits)) {
  sx <- normlower(datTraits$sex)
  datTraits$sex_bin <- ifelse(sx %in% c("male","m"), 1,
                       ifelse(sx %in% c("female","f"), 0, NA_real_))
}

# tissue: colon -> 0, rectum -> 1
if ("tissue" %in% names(datTraits)) {
  tt <- normlower(datTraits$tissue)
  datTraits$tissue_bin <- ifelse(tt %in% c("rectum","rectal"), 1,
                          ifelse(tt %in% c("colon","colonic"), 0, NA_real_))
  attr(datTraits$tissue_bin, "levels") <- c("colon"=0, "rectum"=1)
}

# Keep only numeric columns for WGCNA correlations
is_num <- vapply(datTraits, is.numeric, logical(1))
datTraits_num <- datTraits[, is_num, drop = FALSE]

# Align traits to expression sample order and drop all-NA columns
datTraits_num <- datTraits_num[match(rownames(datExpr0), rownames(datTraits_num)), , drop = FALSE]
all_na <- vapply(datTraits_num, function(v) all(is.na(v)), logical(1))
datTraits_num <- datTraits_num[, !all_na, drop = FALSE]

# Optional: print what traits will be used
message("Traits used for ME correlation: ", paste(colnames(datTraits_num), collapse = ", "))

# Use this for downstream correlations
datTraits <- datTraits_num

#Regrouping samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
#Converting phenotypic characters in a color representation: white means low value, red means high value and gray missing value
traitColors = numbers2colors(datTraits, signed = TRUE)

print("Sample Dendrogram...")
svg("../output/dendrogram_plot.svg", width=8, height=6)
#Plot a sample dendogram with the colors below
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

print("Scale independence and mean connectivity...")
svg("../output/scale_independence_and_mean_connectivity", width=8, height=6)
# Plot scale independence and mean connectivity
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit", type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="blue")
abline(h=0.9, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")
dev.off()

softPower <- 7 # Can change depending on the network
adjacency <- adjacency(datExpr0, power = softPower, type = "signed")

print("Degree distribution")
degree <- rowSums(adjacency)
svg("../output/degree_distribution.svg", width=8, height=6)
hist(degree,
     breaks = 50,
     col = "gray",
     main = "Degree Distribution",
     xlab = "Degree",
     ylab = "Number of Genes")
dev.off()

TOM = TOMsimilarity(adjacency, TOMType = "signed")
dissTOM = 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method="average")

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
                            pamRespectsDendro = FALSE, minClusterSize = 30)

dynamicColors <- labels2colors(dynamicMods)

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs <- MEList$eigengenes

MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

print("Clustering of Module Eigengenes")
svg("../output/clustering_of_module_eigengenes.svg", width=8, height=6)
plot(METree, main="Clustering of Module Eigenes", xlab = "", sub = "")
MEDissThres <- 0.4
abline(h=MEDissThres, col = "red")
dev.off()

merged = mergeCloseModules(datExpr0, dynamicColors, cutHeight = 0.4, verbose = 3)

mergedColors = merged$colors
mergedMEs = merged$newMEs

print("Merged MEs")
svg("../output/merged_MEs.svg", width=8, height=6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors

#Building numeric labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

print("TOMplot")
svg("../output/TOMplot.svg", width=8, height=6)
TOMplot(dissTOM , geneTree,dynamicColors, terrainColors=TRUE)
dev.off()

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

svg("../output/moduletrait_relationships.svg", width=8, height=6)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               colorLabels = FALSE,    # disables color swatches on Y-axis
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

# Convert matrix to a long-format data frame
cor_df <- as.data.frame(as.table(moduleTraitCor))

# Rename columns for clarity
colnames(cor_df) <- c("Module", "Trait", "Correlation")
# Convert p-values matrix to long format
pval_df <- as.data.frame(as.table(moduleTraitPvalue))
colnames(pval_df) <- c("Module", "Trait", "Pvalue")

# Combine correlation and p-value data
cor_df$Pvalue <- pval_df$Pvalue
# Sort by absolute correlation (strongest associations)
cor_df_sorted <- cor_df[order(-abs(cor_df$Correlation)), ]

print(cor_df_sorted)

nSamples <- nrow(datExpr0)
tumoral = as.data.frame(datTraits$PHENOTYPE)
names(tumoral) = "Tumoral"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr0, tumoral, use="p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(tumoral), sep="")
names(GSPvalue) = paste("p.GS", names(tumoral), sep="")

dir.create("../output/mm_vs_gs_plots", showWarnings = FALSE)

# Loop over all modules (skip 'grey' which is unassigned)
for (module in unique(moduleColors)) {
  if (module == "grey") next
  
  column <- match(module, modNames)
  moduleGenes <- moduleColors == module
  
  svg_filename <- paste0("../output/mm_vs_gs_plots/", module, "_membership_vs_gs.svg")
  
  svg(svg_filename, width = 8, height = 6)
  par(mfrow = c(1,1))
  verboseScatterplot(
    abs(geneModuleMembership[moduleGenes, column]),
    abs(geneTraitSignificance[moduleGenes, 1]),
    xlab = paste("Module Membership in", module, "module"),
    ylab = "Gene Significance for Tumoral State",
    main = paste("MM vs GS -", module),
    cex.main = 1.2,
    cex.lab = 1.2,
    cex.axis = 1.2,
    col = module
  )
  dev.off()
}

# make a MDS plot with the expression data, color each gene with the module color
datExpr.genes <- t(datExpr0)

# Compute distance and apply classical MDS
gene.dist <- dist(datExpr.genes)
mds.out <- cmdscale(gene.dist, k = 2)

gene.colors <- labels2colors(moduleColors)

# Basic MDS plot
svg("../output/mds_plot.svg", width=8, height=6)
plot(mds.out,
     col = gene.colors,
     pch = 19,
     main = "MDS plot",
     xlab = "MDS1",
     ylab = "MDS2")
dev.off()

mergedColors <- merged$colors
mergedMEs <- merged$newMEs

# Calculate the dissimilarity of module eigengenes
MEDiss2 <- 1 - cor(mergedMEs)

# Cluster the new module eigengenes
METree2 <- hclust(as.dist(MEDiss2), method = "average")

# Plot the dendrogram
print("Clustering of merged module eigengenes...")
svg("../output/clustering_of_merged_module_eigengenes.svg", width=8, height=6)
plot(METree2, main = "Clustering of Merged Module Eigengenes", xlab = "", sub = "")

# Optional: add the threshold used for merging
abline(h = 0.4, col = "red")
dev.off()

moduleTraitCor <- cor(mergedMEs, datTraits$PHENOTYPE, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr0))
MEs_traits_combined <- cbind(mergedMEs, datTraits$PHENOTYPE)
corMatrix <- cor(MEs_traits_combined, use = "p")

print("Correlation between eigengenes and traits...")
svg("../output/correlation_eigengenes_and_traits.svg", width=8, height=6)
grid.newpage()
p<-pheatmap(corMatrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         main = "Correlation of Module Eigengenes and Traits",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 10,
         fontsize_col = 10)
grid.draw(p$gtable)
dev.off()

write.csv(data.frame(Gene = colnames(datExpr0), Module = moduleColors),
          "../output/gene_module_membership.csv", row.names = FALSE)
module_gene_counts <- as.data.frame(table(moduleColors))
colnames(module_gene_counts) <- c("Module", "GeneCount")
write.csv(module_gene_counts, "../output/module_gene_counts.csv", row.names = FALSE)

degGenes <- rownames(df[df$Category != "Not significant", ])
dir.create("../output/cytoscape_all_modules", showWarnings = FALSE)

for (module in unique(moduleColors)) {
  if (module == "grey") next  # skip unassigned
  
  inModule <- (moduleColors == module)
  modGenes <- colnames(datExpr0)[inModule]
  modGenes_DEG <- intersect(modGenes, degGenes)
  
  # Skip empty modules
  if (length(modGenes_DEG) < 2) next
  
  modExpr <- datExpr0[, modGenes_DEG]
  adjacencyMatrix <- adjacency(modExpr, power = softPower, type = "signed")
  TOM <- TOMsimilarity(adjacencyMatrix)
  rownames(TOM) <- colnames(TOM) <- modGenes_DEG
  
  categoryVec <- df[modGenes_DEG, "Category"]
  categoryVec[is.na(categoryVec)] <- "Not_Annotated"
  
  exportNetworkToCytoscape(
    TOM,
    edgeFile = paste0("../output/cytoscape_all_modules/edges_", module, ".txt"),
    nodeFile = paste0("../output/cytoscape_all_modules/nodes_", module, ".txt"),
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = modGenes_DEG,
    nodeAttr = categoryVec
  )
}

# Here, generate the output that will be input to WGCNA module conservation

dir.create("../output/wgcna_preservation", showWarnings = FALSE)

# Sanity checks
stopifnot(ncol(datExpr0) == length(moduleColors))
stopifnot(identical(colnames(datExpr0), colnames(datExpr0)))  # placeholder to stress gene order matters

# Save one-set export (use one file per cohort/state; you will stitch two of these later)
wgcna_export <- list(
  data      = datExpr0,                 # samples x genes matrix used for WGCNA
  colors    = setNames(moduleColors, colnames(datExpr0)),  # named vector, genes in SAME order
  genes     = colnames(datExpr0),       # explicit gene order manifest
  samples   = rownames(datExpr0),       # sample order manifest
  MEs       = MEs,                      # module eigengenes (optional but handy)
  traits    = datTraits,                # traits aligned to samples (optional)
  softPower = softPower,                # for reference
  networkType = "signed"
)

saveRDS(wgcna_export, file = "../output/wgcna_preservation/export_for_preservation.rds")

message("Saved preservation export to ../output/wgcna_preservation/export_for_preservation.rds")

# ----- Helper (keep in this script for later use) -----
# Given two export RDS files (from two cohorts or two strata),
# this builds multiExpr and colorList aligned on common genes.
build_preservation_inputs <- function(pathA, pathB, refLabel = "SetA", testLabel = "SetB") {
  A <- readRDS(pathA)
  B <- readRDS(pathB)

  # Intersect genes and enforce identical order
  commonGenes <- intersect(A$genes, B$genes)
  if (length(commonGenes) < 1000)
    warning("Few common genes; preservation statistics may be unstable.")

  A.idx <- match(commonGenes, A$genes)
  B.idx <- match(commonGenes, B$genes)

  A.data <- A$data[, A.idx, drop = FALSE]
  B.data <- B$data[, B.idx, drop = FALSE]

  A.cols <- as.integer(A$colors[A.idx])  # colors can be numeric labels or color names; either works
  B.cols <- as.integer(B$colors[B.idx])

  # WGCNA expects a list-of-lists: each set is list(data = samples x genes matrix)
  multiExpr <- list()
  multiExpr[[refLabel]]  <- list(data = A.data)
  multiExpr[[testLabel]] <- list(data = B.data)

  colorList <- list()
  colorList[[refLabel]]  <- A.cols
  colorList[[testLabel]] <- B.cols

  list(multiExpr = multiExpr, colorList = colorList, commonGenes = commonGenes)
}