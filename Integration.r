library(Seurat)
library(SeuratData)
library(patchwork)

# install dataset
InstallData("ifnb")
# This might be hard if the 60s timeout applies

# load dataset
ifnb <- readRDS("C:/Users/ricky/Code/SciP/ifnb.SeuratData/ifnb_processed.Rds")
DefaultAssay(object = ifnb) <- "RNA"

# split the RNA measurements into two layers one for control cells, one for stimulated cells
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb

# run standard analysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)

ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))



# CCA Integration
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", 
                        new.reduction = "integrated.cca")

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap.integrated.cca")
DimPlot(ifnb, reduction = "umap.integrated.cca", group.by = c("stim", "seurat_clusters"))

# Reciprocal PCA
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)

ifnb <- IntegrateLayers(object = ifnb, method = RPCAIntegration, orig.reduction = "pca", 
                        new.reduction = "integrated.rpca")
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.rpca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.integrated.rpca")
DimPlot(ifnb, reduction = "umap.integrated.rpca", group.by = c("stim", "seurat_clusters"))

# Harmony
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)

ifnb <- IntegrateLayers(object = ifnb, method = HarmonyIntegration, orig.reduction = "pca", 
                        new.reduction = "integrated.harmony")
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.harmony", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.harmony", reduction.name = "umap.integrated.harmony")
DimPlot(ifnb, reduction = "umap.integrated.harmony", group.by = c("stim", "seurat_clusters"))
