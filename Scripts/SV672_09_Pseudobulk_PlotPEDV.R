library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(viridis)
library(dplyr)
library(tidyverse)
library(scales)
library(report)
library(edgeR)
library(limma)

# Load Seurat object ----
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples.h5Seurat')
DefaultAssay(seu) <- 'SCT'

# Pseudobulk MDS plot visualization ----
Idents(seu) <- seu$AnimalID
order <- levels(seu$AnimalID)
av.exp <- AverageExpression(seu, return.seurat = FALSE, layer = 'data') # create in-silico bulk RNA-seq dataset for each sample
counts <- as.matrix(av.exp[['SCT']])
colnames(counts) <- gsub('-', '_', colnames(counts))

DGE <- DGEList(counts = counts, genes = rownames(counts), group = colnames(counts)) # make into edgeR DGE object
plotMDS(DGE, 
        top = length(rownames(av.exp[["RNA"]])), #consider all genes 
        col = c('lightblue3', 'steelblue3', 'blue', 'rosybrown2', 'salmon', 'peru', 'red', 'red4'),
        pch = 19,
        cex = 2.5)
legend("bottomright",
       legend = order,
       col = c('lightblue3', 'steelblue3', 'blue', 'rosybrown2', 'salmon', 'peru', 'red', 'red4'),
       pch = 19,
       cex = 1.5)

# PEDV gene expression in each sample
DefaultAssay(seu) <- 'SCT'
VlnPlot(seu, features = 'PEDV', group.by = 'AnimalID', cols = c('lightblue3', 'steelblue3', 'blue', 'rosybrown2', 'salmon', 'peru', 'red', 'red4'))

## Read in counts data ----

seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples_Annotated.h5Seurat')
DefaultAssay(seu) <- 'SCT'
Idents(seu) <- seu$clusters_celltype
seu <- subset(seu, idents = 'B_lowFeature_32', invert = TRUE)
seu$Infection <- paste(seu$Treatment, seu$LocalInfection_IHC, sep = '_')
Idents(seu) <- seu$Infection
seu <- subset(seu, idents = c('Mock_Uninfected', 'PEDV_PEDV'))

## Store percent & scaled expression data ----
seu <- ScaleData(seu, features = rownames(seu))
genes <- c('PEDV')
dat <- as.matrix(seu[["SCT"]]@scale.data[genes,])
seu$combo <- paste(seu$clusters_celltype, seu$Treatment, sep = '_')
dat <- bind_cols(seu$combo, as.data.frame(dat))
colnames(dat) <- c('combo', 'PEDV')
dat <- pivot_longer(dat, -combo, names_to="Gene", values_to="Expression")
dat <- dat %>%
  group_by(combo, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100)
dat$Treatment <- sub('.*\\_', '', dat$combo)
dat$Cluster <- sub("_[^_]+$", "", dat$combo)
dat["Pct"][dat["Pct"] == 0] <- NA
dat$Cluster <- factor(dat$Cluster, levels = rev(c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19', # B_resting
                                                  'B_nonGC_cycling_39', # B_nonGC_cycling
                                                  'B_GC_LZ_2', 'B_GC_LZ_4', # B_GC_LZ
                                                  'B_GC_DZ_1', 'B_GC_DZ_9', # B_GC_DZ
                                                  'ASC_21', 'ASC_27', 'ASC_40', 'ASC_41', # ASC
                                                  'T_ab_resting_12', # T_ab_resting
                                                  'T_CD4_nonnaive_6', # T_CD4_nonnaive
                                                  'T_CD4_follicular_5', # T_CD4_follicular
                                                  'T_CD4_cycling_29', # T_CD4_cycling
                                                  'T_CD8ab_8', 'T_CD8ab_24', # T_CD8ab
                                                  'T_gd_CD2pos_16', 'T_gd_CD2pos_20', # T_gd_CD2pos
                                                  'T_gd_CD2pos_SELLpos_26', # T_gd_CD2pos_SELLpos
                                                  'T_gd_CD2neg_36', # T_gd_CD2neg
                                                  'ILC_group1_ITGAEpos_11', 'ILC_group1_ITGAEpos_23', # ILC_group1_ITGAEpos
                                                  'ILC_group1_ITGAEneg_37', # ILC_group1_ITGAEneg
                                                  'ILC_group3_13', # ILC_group3
                                                  'Macrophage_CD4pos_25', # Macrophage_CD4pos
                                                  'Macrophage_CD4neg_28', # Macrophage_CD4neg
                                                  'cDC_18', # cDC
                                                  'Mast_42', # Mast
                                                  'ISC_TA_22', # ISC_TA
                                                  'Enterocyte_early_17', # Enterocyte_early
                                                  'Enterocyte_intermediate_10', # Enterocyte_intermediate
                                                  'Enterocyte_mature_3', 'Enterocyte_mature_38', # Enterocyte_mature
                                                  'Enterocyte_BEST4pos_33', # Enterocyte_BEST4pos
                                                  'Goblet_30', 'Goblet_34', # Goblet
                                                  'Endothelial_35', # Endothelial
                                                  'Fibroblast_31')))

# Plot genes of interest ----
gene <- genes[1]
#gene <- genes[2]
#gene <- genes[3] # and so on for length(genes)
sub <- subset(dat, Gene == gene)
ggplot(sub, aes(x=Treatment, y = Cluster, color = Avg, size = Pct)) + 
  geom_point() + 
  scale_color_viridis(option = 'rocket', name = 'Average Expression', direction = -1, end = 0.8, begin = 0.2,
                      limits = c(-1, 3), oob = squish) + 
  ylab('') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1, size = 12, face = 'bold'),
        axis.text.y = element_text(size = 12, face = 'bold')) +
  scale_size(range=c(0,6), limits = c(0,100), breaks=c(0,25,50,75,100)) +
  ggtitle(gene)

# Feature plot ----
FeaturePlot(seu, features = 'PEDV', split.by = 'Treatment', reduction = 'umap.cca', order = TRUE, cols = c('grey80', 'darkmagenta'), pt.size = 1)

### View session information ----
report(sessionInfo())
