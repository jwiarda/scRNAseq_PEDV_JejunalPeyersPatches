library(Seurat)
library(SeuratObject)
library(ggplot2)
library(SeuratDisk)
library(SingleCellExperiment)
library(slingshot)
library(viridis)
library(dplyr)
library(tidyverse)
library(tradeSeq)
library(writexl)
library(report)

# Re-process data into subset ----

seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples_Annotated.h5Seurat')
DefaultAssay(seu) <- 'SCT'
seu$Infection <- paste(seu$Treatment, seu$LocalInfection_IHC, sep = '_')
Idents(seu) <- seu$Infection
seu <- subset(seu, idents = c('Mock_Uninfected', 'PEDV_PEDV'))
Idents(seu) <- seu$celllineage
seu <- subset(seu, idents = c('Epithelial'))
DefaultAssay(seu) <- 'SCT'
Idents(seu) <- seu$clusters_celltype

DimPlot(seu, group.by = 'clusters_celltype')

## Re-run dimensionality reductions with just data subset cells ----

seu <- RunPCA(seu, npcs = 50, verbose = F)
ElbowPlot(seu,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... 

seu <- IntegrateLayers(
  object = seu, method = CCAIntegration,
  k.weight = 80, # resolves error: Error: k.weight (100) is set larger than the number of cells in the smallest object (80). Please choose a smaller k.weight.
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = TRUE, normalization.method="SCT"
)

# Slingshot trajectory analysis - all epithelial cells ---

## Convert to a SingleCellExperiment object ----

sc <- as.SingleCellExperiment(seu, assay = 'SCT')
# add original dimensionality reductions to sce object
# could opt to re-calculate these for each data subset if desired
reducedDims(sc) <- SimpleList(CCA=Embeddings(seu@reductions$integrated.cca)[,1:30],
                              UMAP = Embeddings(seu, 'umap.cca'))

#Run Slingshot trajectory analysis:

sc <- slingshot(sc, reducedDim = 'CCA', clusterLabels = 'clusters_celltype', start.clus = 'ISC_TA_22')  # leaving clusters unspecified would allow only a single trajectory to be fit to all cells
summary(sc$slingPseudotime_1) # see summary of pseudotime.... may have more than one trajectory to look at (e.g. $slingPseudotime_2) if clusters were specified in the slingshot() command
summary(sc$slingPseudotime_2)
summary(sc$slingPseudotime_3) 
summary(sc$slingPseudotime_4) # and so on for as many trajectories as are made

## Plot trajectories ----

seu$pseudotime1 <- sc$slingPseudotime_1
seu$pseudotime2 <- sc$slingPseudotime_2
seu$pseudotime3 <- sc$slingPseudotime_3
seu$pseudotime4 <- sc$slingPseudotime_4

### PCA plots ----
seu[['pca.cca']] <- CreateDimReducObject(embeddings = Embeddings(seu@reductions$integrated.cca)[,1:30], key = "PCA_", global = T, assay = "SCT")

notcells <- colnames(seu)[is.na(seu$pseudotime1)]
FeaturePlot(seu, features = c('pseudotime1'), order = TRUE, reduction = 'pca.cca') + 
  scale_colour_viridis(option = 'inferno', begin = 0.1, direction = -1, na.value="grey85") 

notcells <- colnames(seu)[is.na(seu$pseudotime2)]
FeaturePlot(seu, features = c('pseudotime2'), order = TRUE, reduction = 'pca.cca') + 
  scale_colour_viridis(option = 'inferno', begin = 0.1, direction = -1, na.value="grey85") 

notcells <- colnames(seu)[is.na(seu$pseudotime3)]
FeaturePlot(seu, features = c('pseudotime3'), order = TRUE, reduction = 'pca.cca') + 
  scale_colour_viridis(option = 'inferno', begin = 0.1, direction = -1, na.value="grey85") 

notcells <- colnames(seu)[is.na(seu$pseudotime4)]
FeaturePlot(seu, features = c('pseudotime4'), order = TRUE, reduction = 'pca.cca') + 
  scale_colour_viridis(option = 'inferno', begin = 0.1, direction = -1, na.value="grey85") 

DimPlot(seu, reduction = 'pca.cca', group.by = 'clusters_celltype', shuffle = TRUE,
        cols = c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3'))

### Stacked plots ----

dat <- data.frame(seu$pseudotime1, seu$pseudotime2, seu$pseudotime3, seu$pseudotime4, seu$clusters_celltype)
colnames(dat) <- c('pseudotime1', 'pseudotime2', 'pseudotime3', 'pseudotime4', 'cluster')

# trajectory 1:
max <- max(seu$pseudotime1, na.rm = T)
dat <- dat %>% mutate(psuedotime1_bin = cut(pseudotime1, 
                                            breaks=c(0*max, 
                                                     0.1*max,
                                                     0.2*max,
                                                     0.3*max,
                                                     0.4*max,
                                                     0.5*max,
                                                     0.6*max,
                                                     0.7*max,
                                                     0.8*max,
                                                     0.9*max,
                                                     1*max)))
dat$psuedotime1_bin <- as.character(dat$psuedotime1_bin)
dat$psuedotime1_bin <- replace_na(dat$psuedotime1_bin, 'Non-trajectory')

tab <- prop.table(table(dat$cluster, dat$psuedotime1_bin), margin = 2)
tab <- tab[, c('(0,23.2]', '(23.2,46.5]', '(46.5,69.7]', '(69.7,92.9]', '(92.9,116]',
               '(116,139]', '(139,163]', '(163,186]', '(186,209]', '(209,232]')]
tab <- data.frame(tab)
colnames(tab) <- c('Cluster', 'Bin', 'Freq')
tab$Cluster <- factor(tab$Cluster , levels=c('ISC_TA_22', 'Enterocyte_early_17', 'Enterocyte_intermediate_10',
                                             'Enterocyte_mature_3', 'Enterocyte_mature_38', 
                                             'Enterocyte_BEST4pos_33', 'Goblet_30', 'Goblet_34') )
tab <- na.omit(tab)
tab$Bin <- as.numeric(tab$Bin)

ggplot(tab, aes(x=Bin, y=Freq, fill=Cluster)) + geom_area() +
  scale_fill_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

# trajectory 2:
max <- max(seu$pseudotime2, na.rm = T)
dat <- dat %>% mutate(psuedotime2_bin = cut(pseudotime2, 
                                            breaks=c(0*max, 
                                                     0.1*max,
                                                     0.2*max,
                                                     0.3*max,
                                                     0.4*max,
                                                     0.5*max,
                                                     0.6*max,
                                                     0.7*max,
                                                     0.8*max,
                                                     0.9*max,
                                                     1*max)))
dat$psuedotime2_bin <- as.character(dat$psuedotime2_bin)
dat$psuedotime2_bin <- replace_na(dat$psuedotime2_bin, 'Non-trajectory')

tab <- prop.table(table(dat$cluster, dat$psuedotime2_bin), margin = 2)
tab <- tab[, c('(0,19.1]', '(19.1,38.2]', '(38.2,57.3]', '(57.3,76.4]', '(76.4,95.5]',
               '(95.5,115]', '(115,134]', '(134,153]', '(153,172]', '(172,191]')]
tab <- data.frame(tab)
colnames(tab) <- c('Cluster', 'Bin', 'Freq')
tab$Cluster <- factor(tab$Cluster , levels=c('ISC_TA_22', 'Enterocyte_early_17', 'Enterocyte_intermediate_10',
                                             'Enterocyte_mature_3', 'Enterocyte_mature_38', 
                                             'Enterocyte_BEST4pos_33', 'Goblet_30', 'Goblet_34') )
tab <- na.omit(tab)
tab$Bin <- as.numeric(tab$Bin)

ggplot(tab, aes(x=Bin, y=Freq, fill=Cluster)) + geom_area() +
  scale_fill_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

# trajectory 3:
max <- max(seu$pseudotime3, na.rm = T)
dat <- dat %>% mutate(psuedotime3_bin = cut(pseudotime3, 
                                            breaks=c(0*max, 
                                                     0.1*max,
                                                     0.2*max,
                                                     0.3*max,
                                                     0.4*max,
                                                     0.5*max,
                                                     0.6*max,
                                                     0.7*max,
                                                     0.8*max,
                                                     0.9*max,
                                                     1*max)))
dat$psuedotime3_bin <- as.character(dat$psuedotime3_bin)
dat$psuedotime3_bin <- replace_na(dat$psuedotime3_bin, 'Non-trajectory')

tab <- prop.table(table(dat$cluster, dat$psuedotime3_bin), margin = 2)
tab <- tab[, c('(0,39.9]', '(39.9,79.8]', '(79.8,120]', '(120,160]', '(160,199]',
               '(199,239]', '(239,279]', '(279,319]', '(319,359]', '(359,399]')]
tab <- data.frame(tab)
colnames(tab) <- c('Cluster', 'Bin', 'Freq')
tab$Cluster <- factor(tab$Cluster , levels=c('ISC_TA_22', 'Enterocyte_early_17', 'Enterocyte_intermediate_10',
                                             'Enterocyte_mature_3', 'Enterocyte_mature_38', 
                                             'Enterocyte_BEST4pos_33', 'Goblet_30', 'Goblet_34') )
tab <- na.omit(tab)
tab$Bin <- as.numeric(tab$Bin)

ggplot(tab, aes(x=Bin, y=Freq, fill=Cluster)) + geom_area() +
  scale_fill_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

# trajectory 4:
max <- max(seu$pseudotime4, na.rm = T)
dat <- dat %>% mutate(psuedotime4_bin = cut(pseudotime4, 
                                            breaks=c(0*max, 
                                                     0.1*max,
                                                     0.2*max,
                                                     0.3*max,
                                                     0.4*max,
                                                     0.5*max,
                                                     0.6*max,
                                                     0.7*max,
                                                     0.8*max,
                                                     0.9*max,
                                                     1*max)))
dat$psuedotime4_bin <- as.character(dat$psuedotime4_bin)
dat$psuedotime4_bin <- replace_na(dat$psuedotime4_bin, 'Non-trajectory')

tab <- prop.table(table(dat$cluster, dat$psuedotime4_bin), margin = 2)
tab <- tab[, c('(0,28.2]', '(28.2,56.5]', '(56.5,84.7]', '(84.7,113]', '(113,141]',
               '(141,169]', '(169,198]', '(198,226]', '(226,254]', '(254,282]')]
tab <- data.frame(tab)
colnames(tab) <- c('Cluster', 'Bin', 'Freq')
tab$Cluster <- factor(tab$Cluster , levels=c('ISC_TA_22', 'Enterocyte_early_17', 'Enterocyte_intermediate_10',
                                             'Enterocyte_mature_3', 'Enterocyte_mature_38', 
                                             'Enterocyte_BEST4pos_33', 'Goblet_30', 'Goblet_34') )
tab <- na.omit(tab)
tab$Bin <- as.numeric(tab$Bin)

ggplot(tab, aes(x=Bin, y=Freq, fill=Cluster)) + geom_area() +
  scale_fill_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

### Plot PEDV expression ----

FeaturePlot(seu, features = 'PEDV', split.by = 'Treatment', reduction = 'pca.cca',order = TRUE, cols = c('grey80', 'darkmagenta'), pt.size = 1)

dat <- data.frame(seu$pseudotime1, seu$pseudotime2, seu$pseudotime3, seu$pseudotime4, seu$clusters_celltype, seu[["SCT"]]@data['PEDV',], seu$Treatment)
colnames(dat) <- c('pseudotime1', 'pseudotime2', 'pseudotime3', 'pseudotime4', 'cluster', 'PEDV', 'Treatment')
dat <- subset(dat, Treatment == 'PEDV')

sub <- dat[,c('pseudotime1', 'cluster', 'PEDV')]
sub <- na.omit(sub)
ggplot(sub, aes(pseudotime1, PEDV, color = cluster)) + 
  geom_point() +
  scale_color_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() +
  geom_vline(xintercept = biknot, cex = 0.5, color = 'red', linetype='dashed') +
  geom_smooth(aes(color = NULL), color = 'black') # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line

sub <- dat[,c('pseudotime2', 'cluster', 'PEDV')]
sub <- na.omit(sub)
ggplot(sub, aes(pseudotime2, PEDV, color = cluster)) + 
  geom_point() +
  scale_color_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() +
  geom_vline(xintercept = biknot, cex = 0.5, color = 'red', linetype='dashed') +
  geom_smooth(aes(color = NULL), color = 'black') # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line

sub <- dat[,c('pseudotime3', 'cluster', 'PEDV')]
sub <- na.omit(sub)
ggplot(dat, aes(pseudotime3, PEDV, color = cluster)) + 
  geom_point() +
  scale_color_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() +
  geom_smooth(aes(color = NULL), color = 'black') # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line

sub <- dat[,c('pseudotime4', 'cluster', 'PEDV')]
sub <- na.omit(sub)
ggplot(dat, aes(pseudotime4, PEDV, color = cluster)) + 
  geom_point() +
  scale_color_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() +
  geom_smooth(aes(color = NULL), color = 'black') # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line


# Save files ----
dir.create('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/')
saveRDS(sc, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/slingshot_Epithelial.rds')
SaveH5Seurat(seu, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/trajectorySubset_Epithelial.h5seurat')

# Define knots ----

# first identify number of knots to use
genes <- VariableFeatures(seu) # find most variable genes to test...not testing all genes will severely cut down on run time
set.seed(123) # try this with a couple seeds to make sure results were reproducible
icMat <- evaluateK(counts = as.matrix(assays(sc)$counts[genes,]),
                   pseudotime = as.matrix(sc$slingshot@assays@data$pseudotime),
                   cellWeights = as.matrix(sc$slingshot@assays@data$weights),
                   conditions = factor(colData(sc)$Treatment),
                   nGenes = 300,
                   k = 3:10)
nknots = 7 # use 7 knots for this data; tend to choose lower number as to not overfit; choose knot with little variation or where variation begins to level off in middle plots

# Set GAM model ----

set.seed(123)
BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 4 
sc <- fitGAM(counts = as.matrix(assays(sc)$counts[genes,]),
             sds = sc$slingshot,
             conditions = factor(colData(sc)$Treatment),
             nknots = nknots,
             BPPARAM = BPPARAM,
             parallel = TRUE,
             verbose = TRUE)
mean(rowData(sc)$tradeSeq$converged)
biknot <- sc@metadata$tradeSeq$knots[4]
saveRDS(sc, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/pseudotime_GAM_model_Epithelial.rds')

# Subset trajectory - polarized enterocytes ----
#sc@metadata$tradeSeq$knots[4] <- 148.7960
dat <- data.frame(seu$pseudotime1, seu$pseudotime2, seu$clusters_celltype)
dat <- dat[!with(dat,is.na(seu.pseudotime1) & is.na(seu.pseudotime2)),]
dat[is.na(dat)] <- -5

# identify knot value of bifurcation ----
ggplot(dat, aes(seu.pseudotime1, seu.pseudotime2, color = seu.clusters_celltype)) + 
  geom_point(size = 0.2, position = position_jitter(width = 5, height = 5, seed = 100)) +
  scale_color_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() +
  geom_vline(xintercept = biknot, cex = 0.5, color = 'red', linetype='dashed') +
  geom_hline(yintercept = biknot, cex = 0.5, color = 'red', linetype='dashed') +
  geom_point(aes(x=biknot, y=biknot), colour="black", size = 4, shape = 8)

t1 <- subset(dat, seu.pseudotime1 > biknot)
t2 <- subset(dat, seu.pseudotime2 > biknot)
seu3 <- subset(seu, cells = c(rownames(t1), rownames(t2)))

DefaultAssay(seu3) <- 'SCT'
Idents(seu3) <- seu3$clusters_celltype

## Re-run dimensionality reductions with just data subset cells ----

seu3 <- RunPCA(seu3, npcs = 50, verbose = F)
ElbowPlot(seu3,
          ndims = 50) # look at this plot to find the 'elbow' for significant PCs... 

seu3 <- IntegrateLayers(
  object = seu3, method = CCAIntegration, 
  k.weight = 36, # resolves error: Error: k.weight (100) is set larger than the number of cells in the smallest object (36). Please choose a smaller k.weight.
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = TRUE, normalization.method="SCT"
)

# Slingshot trajectory analysis ---

## Convert to a SingleCellExperiment object ----

sc <- as.SingleCellExperiment(seu3, assay = 'SCT')
# add original dimensionality reductions to sce object
# could opt to re-calculate these for each data subset if desired
reducedDims(sc) <- SimpleList(CCA=Embeddings(seu3@reductions$integrated.cca)[,1:30],
                              UMAP = Embeddings(seu3, 'umap.cca'))

#Run Slingshot trajectory analysis:

sc <- slingshot(sc, reducedDim = 'CCA')  # leaving clusters unspecified would allow only a single trajectory to be fit to all cells
summary(sc$slingPseudotime_1) # see summary of pseudotime.... may have more than one trajectory to look at (e.g. $slingPseudotime_2) if clusters were specified in the slingshot() command
#summary(sc$slingPseudotime_2) # and so on for as many trajectories as are made

## Plot trajectories ----

seu3$pseudotime1new <- sc$slingPseudotime_1

# since a start vs end of trajectory wasn't established, we could view either end as start. I prefer to look at Traj 1 cells as start (though vice versa could be true) so inverted trajectory pseudotime values
max(seu3$pseudotime1new)
min(seu3$pseudotime1new)
seu3$pseudotime1new <- (seu3$pseudotime1new-max(seu3$pseudotime1new))*(-1)

### PCA plots ----
seu3[['pca.cca']] <- CreateDimReducObject(embeddings = Embeddings(seu3@reductions$integrated.cca)[,1:30], key = "PCA_", global = T, assay = "SCT")

notcells <- colnames(seu3)[is.na(seu3$pseudotime1new)]
#FeaturePlot(seu3, features = c('pseudotime1new'), order = TRUE, reduction = 'pca.cca') + 
#  scale_colour_viridis(option = 'inferno', begin = 0.1, direction = -1, na.value="grey85") 
FeaturePlot(seu3, features = c('pseudotime1new'), order = TRUE, reduction = 'pca.cca', cols = c('limegreen', 'purple1'))

notcells <- colnames(seu3)[is.na(seu3$pseudotime1)]
FeaturePlot(seu3, features = c('pseudotime1'), order = TRUE, reduction = 'pca.cca') + 
  scale_colour_viridis(option = 'inferno', begin = 0.1, direction = -1, na.value="grey85") 

notcells <- colnames(seu3)[is.na(seu3$pseudotime2)]
FeaturePlot(seu3, features = c('pseudotime2'), order = TRUE, reduction = 'pca.cca') + 
  scale_colour_viridis(option = 'inferno', begin = 0.1, direction = -1, na.value="grey85") 

DimPlot(seu3, reduction = 'pca.cca', group.by = 'clusters_celltype', shuffle = TRUE,
        cols = c('cyan4', 'gold3', 'chartreuse4', 'darkmagenta'))

### Stacked plots ----

# trajectory 1 cell types:
dat <- data.frame(seu3$pseudotime1, seu3$pseudotime2, seu3$pseudotime1new, seu3$clusters_celltype)
colnames(dat) <- c('pseudotime1', 'pseudotime2', 'pseudotime1new', 'cluster')
max <- max(seu3$pseudotime1new, na.rm = T)
dat <- dat %>% mutate(psuedotime1_bin = cut(pseudotime1new, 
                                            breaks=c(0*max, 
                                                     #0.1*max,
                                                     0.2*max,
                                                     #0.3*max,
                                                     0.4*max,
                                                     #0.5*max,
                                                     0.6*max,
                                                     #0.7*max,
                                                     0.8*max,
                                                     #0.9*max,
                                                     1*max)))
dat$psuedotime1_bin <- as.character(dat$psuedotime1_bin)
dat$psuedotime1_bin <- replace_na(dat$psuedotime1_bin, 'Non-trajectory')

tab <- prop.table(table(dat$cluster, dat$psuedotime1_bin), margin = 2)
tab <- tab[, c('(0,33.7]', '(33.7,67.4]', '(67.4,101]', '(101,135]', '(135,168]')]
tab <- data.frame(tab)
colnames(tab) <- c('Cluster', 'Bin', 'Freq')
tab$Cluster <- factor(tab$Cluster , levels=c('ISC_TA_22', 'Enterocyte_early_17', 'Enterocyte_intermediate_10',
                                             'Enterocyte_mature_3', 'Enterocyte_mature_38', 
                                             'Enterocyte_BEST4pos_33', 'Goblet_30', 'Goblet_34') )
tab <- na.omit(tab)
tab$Bin <- as.numeric(tab$Bin)

ggplot(tab, aes(x=Bin, y=Freq, fill=Cluster)) + geom_area() +
  scale_fill_manual(values=c('sandybrown', 'cyan4', 'gold3', 'chartreuse4', 'darkmagenta', 'hotpink', 'firebrick', 'lightskyblue3')) +
  theme_bw() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

# trajectory 1 polarized states:
dat <- data.frame(seu3$pseudotime1, seu3$pseudotime2, seu3$pseudotime1new, seu3$clusters_celltype)
colnames(dat) <- c('pseudotime1', 'pseudotime2', 'pseudotime1new', 'cluster')
dat[is.na(dat)] <- -5
max <- max(seu3$pseudotime1new, na.rm = T)
dat <- dat %>% mutate(psuedotime1_bin = cut(pseudotime1new, 
                                            breaks=c(0*max, 
                                                     #0.1*max,
                                                     0.2*max,
                                                     #0.3*max,
                                                     0.4*max,
                                                     #0.5*max,
                                                     0.6*max,
                                                     #0.7*max,
                                                     0.8*max,
                                                     #0.9*max,
                                                     1*max)))
dat$psuedotime1_bin <- as.character(dat$psuedotime1_bin)
dat$psuedotime1_bin <- replace_na(dat$psuedotime1_bin, 'Non-trajectory')

dat$polar <- 'Both'
dat$polar[dat$pseudotime1 > 0 & dat$pseudotime2 < 0] <- 'Traj1'
dat$polar[dat$pseudotime2 > 0 & dat$pseudotime1 < 0] <- 'Traj2'

tab <- prop.table(table(dat$polar, dat$psuedotime1_bin), margin = 2)
tab <- tab[, c('(0,33.7]', '(33.7,67.4]', '(67.4,101]', '(101,135]', '(135,168]')]
tab <- data.frame(tab)
colnames(tab) <- c('Cluster', 'Bin', 'Freq')
tab$Cluster <- factor(tab$Cluster , levels=c('Traj1', 'Both', 'Traj2') )
tab <- na.omit(tab)
tab$Bin <- as.numeric(tab$Bin)

ggplot(tab, aes(x=Bin, y=Freq, fill=Cluster)) + geom_area() +
  scale_fill_manual(values=c('limegreen', 'grey50', 'purple1')) +
  theme_bw() + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

### Plot PEDV expression ----

FeaturePlot(seu3, features = 'PEDV', split.by = 'Treatment', reduction = 'pca.cca',order = TRUE, cols = c('grey80', 'darkmagenta'), pt.size = 1)

dat <- data.frame(seu3$pseudotime1, seu3$pseudotime2, seu3$pseudotime1new, seu3$clusters_celltype, seu3[["SCT"]]@data['PEDV',], seu3$Treatment)
colnames(dat) <- c('pseudotime1', 'pseudotime2', 'pseudotime1new', 'cluster', 'PEDV', 'Treatment')
dat <- subset(dat, Treatment == 'PEDV')

sub <- dat[,c('pseudotime1new', 'cluster', 'PEDV')]
sub <- na.omit(sub)
ggplot(sub, aes(pseudotime1new, PEDV, color = cluster)) + 
  geom_point() +
  scale_color_manual(values=c('cyan4', 'gold3', 'chartreuse4', 'darkmagenta')) +
  theme_bw() +
  geom_smooth(aes(color = NULL), color = 'black') # ADD A LOESS SMOOTH LINE with default span of 0.75 and confidence interval shown around line

# Save files ----

#dir.create('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/')
saveRDS(sc, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/slingshot_PolarizedEnterocytes.rds')
SaveH5Seurat(seu3, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/trajectorySubset_PolarizedEnterocytes.h5seurat')

### View session information ----
report(sessionInfo())
