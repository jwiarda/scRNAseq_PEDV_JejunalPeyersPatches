library(Seurat)
library(future)
library(SeuratDisk)
library(clustree)
library(ggplot2)
library(dplyr)
library(AUCell)
library(plyr)
library(report)

seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples_Annotated.h5Seurat')
DefaultAssay(seu) <- 'SCT'

Idents(seu) <- seu$clusters_celltype
seu <- subset(seu, idents = 'B_lowFeature_32', invert = TRUE)
seu$Infection <- paste(seu$Treatment, seu$LocalInfection_IHC, sep = '_')
Idents(seu) <- seu$Infection
seu <- subset(seu, idents = c('Mock_Uninfected', 'PEDV_PEDV'))

seu$combo <- paste(seu$clusters_celltype, seu$Treatment, sep = '_')

# Identify gene signature ----
## Read in & filter DE data ----
AllCellsDGE <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_ClustersPairwise_PEDVEnterocyteInfectionStates_unflitered.tsv')
AllCellsDGE <- subset(AllCellsDGE, logFC > 0.25 & padj < 0.05)
AllCellsDGE <- subset(AllCellsDGE, pct_in > 0.1 | pct_out > 0.1)

Hs_BystanderVsMock <- subset(AllCellsDGE, group == 'HomeostaticEnterocyte_Bystander' & (pop1 == 'HomeostaticEnterocyte_Mock' | pop2 == 'HomeostaticEnterocyte_Mock'))
Hs_BystanderVsMock$group <- 'Hs_BystanderVsMock'
Hs_InfectedVsMock <- subset(AllCellsDGE, group == 'HomeostaticEnterocyte_Infected' & (pop1 == 'HomeostaticEnterocyte_Mock' | pop2 == 'HomeostaticEnterocyte_Mock'))
Hs_InfectedVsMock$group <- 'Hs_InfectedVsMock'

St_BystanderVsMock <- subset(AllCellsDGE, group == 'StressedEnterocyte_Bystander' & (pop1 == 'StressedEnterocyte_Mock' | pop2 == 'StressedEnterocyte_Mock'))
St_BystanderVsMock$group <- 'St_BystanderVsMock'
St_InfectedVsMock <- subset(AllCellsDGE, group == 'StressedEnterocyte_Infected' & (pop1 == 'StressedEnterocyte_Mock' | pop2 == 'StressedEnterocyte_Mock'))
St_InfectedVsMock$group <- 'St_InfectedVsMock'

ConservedPEDVSignature <- intersect(St_InfectedVsMock$feature, (intersect(St_BystanderVsMock$feature, (intersect(Hs_BystanderVsMock$feature, Hs_InfectedVsMock$feature)))))
length(ConservedPEDVSignature)

# Perform gene set enrichment analysis ----
geneSets <- list(InfectionSignature = ConservedPEDVSignature)
exprMatrix <- as.matrix(seu[["SCT"]]@counts) 
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) # calculate cell signature AUC score for each gene set in each cell

clusterID <- as.data.frame(seu[["combo"]])
clusterID <- t(clusterID)
AUCs <- as.data.frame(getAUC(cells_AUC))
AUCs <- rbind(AUCs, clusterID)
AUCs <- t(AUCs)
AUCs <- as.data.frame(AUCs)
head(AUCs)

# Plot AUC scores----
geneSetName <- rownames(cells_AUC)[grep("InfectionSignature", rownames(cells_AUC))]
thres = 0.35
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=thres)

AUCs$InfectionSignature <- as.numeric(AUCs$InfectionSignature)
seu$InfectionSignature <- AUCs$InfectionSignature

FeaturePlot(seu, features = 'InfectionSignature', reduction = 'umap.cca', cols = c('gold', 'navy'), split.by = 'Treatment')

VlnPlot(seu, features = 'InfectionSignature', pt.size = 0.00, group.by = 'clusters_celltype', split.by = 'Treatment', cols = c('dodgerblue3', 'red3'))

# Plot against trajectories ----
## Epithelial differentiation trajectories
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/trajectorySubset_Epithelial.h5seurat')
DefaultAssay(seu) <- 'SCT'

AUCs <- AUCs[rownames(AUCs) %in% colnames(seu),]
identical(rownames(AUCs), colnames(seu)) # make sure TRUE
seu <- AddMetaData(seu, AUCs)

dat <- data.frame(seu$pseudotime1, seu$pseudotime2, seu$pseudotime3, seu$pseudotime4, seu$clusters_celltype, seu$InfectionSignature, seu$Treatment)
colnames(dat) <- c('pseudotime1', 'pseudotime2', 'pseudotime3', 'pseudotime4', 'cluster', 'InfectionSignature', 'Treatment')
#dat <- subset(dat, Treatment == 'PEDV')

sub <- dat[,c('pseudotime1', 'cluster', 'InfectionSignature', 'Treatment')]
sub <- na.omit(sub)
ggplot(sub, aes(pseudotime1, InfectionSignature)) + 
  geom_point(aes(color = cluster)) +
  scale_color_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() +
  geom_smooth(aes(group = Treatment), data = filter(sub, Treatment == 'PEDV'), color = 'red', linewidth = 2) + # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line
  geom_smooth(aes(group = Treatment), data = filter(sub, Treatment == 'Mock'), color = 'blue', linewidth = 2)  # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line

sub <- dat[,c('pseudotime2', 'cluster', 'InfectionSignature', 'Treatment')]
sub <- na.omit(sub)
ggplot(sub, aes(pseudotime2, InfectionSignature)) + 
  geom_point(aes(color = cluster)) +
  scale_color_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() +
  geom_smooth(aes(group = Treatment), data = filter(sub, Treatment == 'PEDV'), color = 'red', linewidth = 2) + # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line
  geom_smooth(aes(group = Treatment), data = filter(sub, Treatment == 'Mock'), color = 'blue', linewidth = 2)  # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line

## Polarized enterocyte transition trajectory
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/trajectorySubset_PolarizedEnterocytes.h5seurat')
DefaultAssay(seu) <- 'SCT'

AUCs <- AUCs[rownames(AUCs) %in% colnames(seu),]
identical(rownames(AUCs), colnames(seu)) # make sure TRUE
seu <- AddMetaData(seu, AUCs)

dat <- data.frame(seu$pseudotime1new, seu$clusters_celltype, seu$InfectionSignature, seu$Treatment)
colnames(dat) <- c('pseudotime1new', 'cluster', 'InfectionSignature', 'Treatment')
#dat <- subset(dat, Treatment == 'PEDV')

sub <- dat[,c('pseudotime1new', 'cluster', 'InfectionSignature', 'Treatment')]
sub <- na.omit(sub)
ggplot(sub, aes(pseudotime1new, InfectionSignature)) + 
  geom_point(aes(color = cluster)) +
  scale_color_manual(values=c('cyan4', 'gold3', 'chartreuse4', 'darkmagenta')) +
  theme_bw() +
  geom_smooth(aes(group = Treatment), data = filter(sub, Treatment == 'PEDV'), color = 'red3') + # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line
  geom_smooth(aes(group = Treatment), data = filter(sub, Treatment == 'Mock'), color = 'dodgerblue3')  # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line

### View session information ----
report(sessionInfo())
