library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(scales)
library(report)
library(glmGamPoi) # speeds up SCTranform
library(future)
library(harmony)
plan("multisession", workers = 4)
options(future.globals.maxSize = 15000 * 1024^2)

# Load Seurat object ----
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/FilterQCOutputs/AllSamples_Filtered.h5Seurat')

#See how many cells/genes in Seurat object:
seu

#See how many cells in each sample:
table(seu$orig.ident)

# SCTransform ----

#options(future.globals.maxSize = 3e+09)
#options(future.globals.maxSize = 10000 * 1024^2)
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident) # variable to use for integration
seu <- SCTransform(seu)
seu # check SCT is present and the active assay


# PCA ----

seu <- RunPCA(seu, npcs = 50, verbose = F)

#Visualize PCs:
ElbowPlot(seu,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... 

plan("multisession", workers = 6)
options(future.globals.maxSize = 30000 * 1024^2)

# Integrate ----

seu <- IntegrateLayers(
  object = seu, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = TRUE, normalization.method="SCT"
)

seu <- IntegrateLayers(
  object = seu, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = TRUE, normalization.method="SCT"
)

seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE, normalization.method="SCT"
)

# Dim reduction ----

seu <- RunUMAP(seu, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
seu <- RunTSNE(seu, reduction = "integrated.cca", dims = 1:30, reduction.name = "tsne.cca")
seu <- RunUMAP(seu, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
seu <- RunTSNE(seu, reduction = "integrated.rpca", dims = 1:30, reduction.name = "tsne.rpca")
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
seu <- RunTSNE(seu, reduction = "harmony", dims = 1:30, reduction.name = "tsne.harmony")
seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, reduction.name = "umap.pca_sct")
seu <- RunTSNE(seu, reduction = "pca", dims = 1:30, reduction.name = "tsne.pca_sct")

DimPlot(seu, shuffle=TRUE, group.by = 'orig.ident', reduction = 'umap.cca')
DimPlot(seu, shuffle=TRUE, group.by = 'orig.ident', reduction = 'umap.rpca')
DimPlot(seu, shuffle=TRUE, group.by = 'orig.ident', reduction = 'umap.harmony')
DimPlot(seu, shuffle=TRUE, group.by = 'orig.ident', reduction = 'umap.pca_sct')

# Normalize/scale RNA assay ----

seu <- NormalizeData(seu,  # normalize the RNA counts data per cell
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000, 
                     assay = "RNA")
seu <- ScaleData(seu, # scale the RNA counts data relative to other cells
                 assay = "RNA")

# Save object ----

seu[["RNA"]] <- as(object = seu[["RNA"]], Class = "Assay")
dir.create('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/')
SaveH5Seurat(seu, 
             filename = "/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples.h5Seurat", overwrite = TRUE)

### View session information ----
report(sessionInfo())
