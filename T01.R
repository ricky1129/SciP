#---
# title: "T01"
# output: html_document
#---

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/ricky/Code/SciP/filtered_gene_bc_matrices/hg19") # nolint
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

#Let's see what pbmc is
pbmc


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Select resonable cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.

# We can also normalize by this, which is exactly the same as the parametres are default values. 
pbmc <- NormalizeData(pbmc)



# Feature Selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#Scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)




# Linear dimentional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Plot 2 expressions separately
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# Plot 2 expressions on 2 axes
DimPlot(pbmc, reduction = "pca") + NoLegend()

# Heat Map
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# Multiple Genes
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# Elbow Plot (Dimontionality)
ElbowPlot(pbmc)


# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


# Non-linear dimensional reduction
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")



# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
