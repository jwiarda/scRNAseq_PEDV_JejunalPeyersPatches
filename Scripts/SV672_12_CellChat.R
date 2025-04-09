library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(report)

seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples_Annotated.h5Seurat')
DefaultAssay(seu) <- 'SCT'
Idents(seu) <- seu$clusters_celltype
seu <- subset(seu, idents = 'B_lowFeature_32', invert = TRUE)
seu$Infection <- paste(seu$Treatment, seu$LocalInfection_IHC, sep = '_')
Idents(seu) <- seu$Infection
seu <- subset(seu, idents = c('Mock_Uninfected', 'PEDV_PEDV'))

# Remove cells belonging to cell types not well represented in both treatments: ----
df <- as.data.frame.matrix(table(seu$clusters_celltype, seu$Infection)) # create data frame of cell type quantities per tissue
order <- c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19', # B_resting
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
           'Fibroblast_31') # define order of cell type levels ot use
df <- df[match(order, rownames(df)), ]  # change to order we would like to present annotations in
min(df) # if minimum is 10 or greater, the next steps won't filter out any cell populations
rm <- rownames(df[df$Mock_Uninfected < 10 | df$PEDV_PEDV < 10,]) # identify cell types that don't have at least 10 cells recovered in each tissue & remove these from further analysis
df <- df[!(row.names(df) %in% rm),]
group.new <- rownames(df) # this is a vector of cell types 
Idents(seu) <- seu$clusters_celltype
seu <- subset(seu, idents = c(group.new)) # subset to only cell types with >= 10 cells in each tissue
Idents(seu) <- seu$clusters_celltype
levels(seu) <- order
seu$clusters_celltype <- Idents(seu)
DefaultAssay(seu) <- 'SCT'

# Convert to human ortholog genes: ----
# Load in ortho genes
orthoGenes <- read.delim("/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/PigToHuman_GeneOrthos_v97.txt") # read in gene ortholog file
orthoGenes <- subset(orthoGenes, Human.homology.type == 'ortholog_one2one') # subset to only one to one orthologs

# Filter  data to include only one-to-one gene orthologs & convert to human gene symbols:
genes <- as.data.frame(rownames(seu[['SCT']]@data)) # extract pig gene names from dataset
colnames(genes) <- 'gene'
pigGenes <- read_excel('/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx') # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalList <-gsub("_", "-", pigGenes$FinalList) # replace all underscores with dashes since this occurred when processing data in a previous step
pigGenes <- pigGenes[pigGenes$FinalList %in% genes$gene, ] # slim down to only genes in our dataset
orthos <- intersect(pigGenes$ENSID, orthoGenes$Gene.stable.ID) # find which genes are one-to-one orthologs
#length(orthos) # how many genes are orthologs?
pigGenes <- pigGenes[pigGenes$ENSID %in% orthos, ]
#dim(pigGenes)
orthoGenes <- orthoGenes[orthoGenes$Gene.stable.ID %in% pigGenes$ENSID, ] # slim down to only ortho genes in our dataset
orthoGenes <- orthoGenes %>% distinct(orthoGenes$Gene.stable.ID, orthoGenes$Human.gene.stable.ID, .keep_all = TRUE)  # retain only unique combinations of pig & human Ensembl IDs, ignoring transcript IDs
#dim(orthoGenes) # should not have same number of rows as in pigGenes
counts <- seu[['SCT']]@data[rownames(seu[['SCT']]@data) %in% pigGenes$FinalList,]
pigGenes <- pigGenes %>% arrange(factor(FinalList, levels = rownames(counts))) # arrange pigGenes in same order as counts
orthoGenes <- orthoGenes %>% arrange(factor(Gene.stable.ID, levels = pigGenes$ENSID)) # arrange orthoGenes in same order as pigGenes (and consequently counts)
rownames(counts) <- orthoGenes$Human.gene.name # change pig genes to human gene names

# Do a custom bit here to incorporate more MHC-II SLA genes, since some are one-to-many orthologs but we are highly interested in MHC-II related pathways for this work...
add <- setdiff(rownames(seu[['SCT']]@data)[which(grepl("SLA-D", rownames(seu[['SCT']]@data)))],
               pigGenes$FinalList[which(grepl("SLA-D", pigGenes$FinalList))]) # ID MHC-II genes from Seurat object that weren't recognized as one-to-one orthologs
counts2 <- seu[['SCT']]@data[rownames(seu[['SCT']]@data) %in% add,]
rownames(counts2) <- gsub('SLA-','HLA-', rownames(counts2)) # change to human gene symbol IDs
counts <- rbind(counts, counts2)

# Make data subsets ----
# Subset to treatments separately:
meta <- data.frame(seu@meta.data)
cellsMock <- rownames(subset(meta, Infection == 'Mock_Uninfected'))
cellsPEDV <- rownames(subset(meta, Infection == 'PEDV_PEDV'))

countsMock <- as.matrix(counts[, cellsMock])
countsPEDV <- as.matrix(counts[, cellsPEDV])

meta$samples <-  meta$AnimalID
metaMock <- subset(meta, Infection == 'Mock_Uninfected')
metaPEDV <- subset(meta, Infection == 'PEDV_PEDV')

# All signaling pathways ----

## Create CellChat object: ----
cellchatMock <- createCellChat(object = countsMock, # use new Seurat object with human gene names
                               meta = metaMock,
                               group.by = "clusters_celltype") # set cell identities to cell type annotations
cellchatPEDV <- createCellChat(object = countsPEDV, # use new Seurat object with human gene names
                               meta = metaPEDV,
                               group.by = "clusters_celltype") # set cell identities to cell type annotations

# reset order:
cellchatMock@idents = factor(cellchatMock@idents,
                             levels = order)
cellchatPEDV@idents = factor(cellchatPEDV@idents,
                             levels = order)

## Set cell interaction database: ----
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

CellChatDB.use <- CellChatDB # look at all signaling categories
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling only
#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # look at only cell-cell contact
#CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # only look at extracellular matrix signaling
#CellChatDB.use <- subsetDB(CellChatDB, search = "Non-protein Signaling") # only look at non-protein signals

# set the used database in the object
cellchatMock@DB <- CellChatDB.use
cellchatPEDV@DB <- CellChatDB.use

## Preprocess expression data ----
# subset the expression data of signaling genes for saving computation cost
cellchatMock <- subsetData(cellchatMock) # This step is necessary even if using the whole database
cellchatPEDV <- subsetData(cellchatPEDV) # This step is necessary even if using the whole database

cellchatMock <- identifyOverExpressedGenes(cellchatMock)
cellchatMock <- identifyOverExpressedInteractions(cellchatMock)
cellchatPEDV <- identifyOverExpressedGenes(cellchatPEDV)
cellchatPEDV <- identifyOverExpressedInteractions(cellchatPEDV)

## Compute communication probabilities and network inferences: ----
#cellchat <- computeCommunProb(cellchat)
cellchatMock <- computeCommunProb(cellchatMock, type =  "truncatedMean", trim = 0.1) # count as zero expression if in <10% of annotated cell type (default = 25%)
cellchatPEDV <- computeCommunProb(cellchatPEDV, type =  "truncatedMean", trim = 0.1) # count as zero expression if in <10% of annotated cell type (default = 25%)

## Infer cell signaling pathway communication & calculate aggregated communication network: ----
cellchatMock <- computeCommunProbPathway(cellchatMock)
cellchatMock <- aggregateNet(cellchatMock)
cellchatPEDV <- computeCommunProbPathway(cellchatPEDV)
cellchatPEDV <- aggregateNet(cellchatPEDV)

## Compute network centrality
cellchatMock <- netAnalysis_computeCentrality(cellchatMock, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchatPEDV <- netAnalysis_computeCentrality(cellchatPEDV, slot.name = "netP")

## Save CellChat object: ----
dir.create('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellChat')
saveRDS(cellchatMock, file = '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellChat/cellchatMock_AllSignaling.rds')
saveRDS(cellchatPEDV, file = '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellChat/cellchatPEDV_AllSignaling.rds')


# Read in all datasets ----
rm(list = ls())
cellchatMock_All <- readRDS('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellChat/cellchatMock_AllSignaling.rds')
cellchatPEDV_All <- readRDS('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellChat/cellchatPEDV_AllSignaling.rds')

# Create & save merged datasets ----

## update objects ----
cellchatMock_All <- updateCellChat(cellchatMock_All)
cellchatPEDV_All <- updateCellChat(cellchatPEDV_All)

## merge objects ----
object.list <- list(Mock = cellchatMock_All, PEDV = cellchatPEDV_All)
cellchat_All <- mergeCellChat(object.list, add.names = names(object.list))

## save merged objects ----
saveRDS(cellchat_All, file = '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellChat/cellchat_TreatmentsMerged_AllSignaling.rds')

# Plot individual dataset heatmaps of interaction number/strength ----
par(mfrow=c(1,1))
groupSize <- as.numeric(table(cellchatMock_All@idents))
g1 <- netVisual_heatmap(cellchatMock_All, title.name = 'All Signaling (Strength) - Mock', measure = 'weight', color.heatmap = c('khaki', 'purple4'),
                        color.use = c(rep('gold3', 5), 'deeppink4', rep('chartreuse4', 2),
                                      rep('cyan4', 2), rep('sandybrown', 4), 'navy', 'lightpink', 
                                      'salmon', 'deepskyblue', rep('tan4', 2), rep('mediumpurple1', 2),
                                      'darkgreen', 'gray50', rep('darkmagenta', 2), 'red',
                                      'hotpink', 'khaki', 'limegreen', 'cadetblue', 'firebrick',
                                      'deepskyblue4', 'darkseagreen', 'burlywood3', 
                                      rep('goldenrod3', 2), 'blue', rep('deeppink', 2),
                                      'lightskyblue3', 'mistyrose3'))
groupSize <- as.numeric(table(cellchatPEDV_All@idents))
g2 <- netVisual_heatmap(cellchatPEDV_All, title.name = 'All Signaling (Strength) - PEDV', measure = 'weight', color.heatmap = c('khaki', 'purple4'),
                        color.use = c(rep('gold3', 5), 'deeppink4', rep('chartreuse4', 2),
                                      rep('cyan4', 2), rep('sandybrown', 4), 'navy', 'lightpink', 
                                      'salmon', 'deepskyblue', rep('tan4', 2), rep('mediumpurple1', 2),
                                      'darkgreen', 'gray50', rep('darkmagenta', 2), 'red',
                                      'hotpink', 'khaki', 'limegreen', 'cadetblue', 'firebrick',
                                      'deepskyblue4', 'darkseagreen', 'burlywood3', 
                                      rep('goldenrod3', 2), 'blue', rep('deeppink', 2),
                                      'lightskyblue3', 'mistyrose3'))
g1 + g2


# Heatmap of differential interactions across treatments ----

# red = up in PEDV, blue = up in mock
netVisual_heatmap(cellchat_All, title.name = 'All Signaling Interactions (Strength)', measure = "weight", color.use = c(rep('gold3', 5), 'deeppink4', rep('chartreuse4', 2),
                                                                                                                               rep('cyan4', 2), rep('sandybrown', 4), 'navy', 'lightpink', 
                                                                                                                               'salmon', 'deepskyblue', rep('tan4', 2), rep('mediumpurple1', 2),
                                                                                                                               'darkgreen', 'gray50', rep('darkmagenta', 2), 'red',
                                                                                                                               'hotpink', 'khaki', 'limegreen', 'cadetblue', 'firebrick',
                                                                                                                               'deepskyblue4', 'darkseagreen', 'burlywood3', 
                                                                                                                               rep('goldenrod3', 2), 'blue', rep('deeppink', 2),
                                                                                                                               'lightskyblue3', 'mistyrose3')) # strength of interactions

# Scatter plot of differential signaling across treatments ----

# positive values = PEDV, negative values = mock

netAnalysis_diff_signalingRole_scatter(cellchat_All, title = 'All Signaling Interactions (Strength)', x.measure = "outdeg", y.measure = "indeg", xlabel = "Outgoing interaction strength", ylabel = "Incoming interaction strength",
                                       color.use = c(rep('gold3', 5), 'deeppink4', rep('chartreuse4', 2),
                                                     rep('cyan4', 2), rep('sandybrown', 4), 'navy', 'lightpink', 
                                                     'salmon', 'deepskyblue', rep('tan4', 2), rep('mediumpurple1', 2),
                                                     'darkgreen', 'gray50', rep('darkmagenta', 2), 'red',
                                                     'hotpink', 'khaki', 'limegreen', 'cadetblue', 'firebrick',
                                                     'deepskyblue4', 'darkseagreen', 'burlywood3', 
                                                     rep('goldenrod3', 2), 'blue', rep('deeppink', 2),
                                                     'lightskyblue3', 'mistyrose3'))

# Identify pathways with largest difference in information flows between treatments ----
## APCs --> Tfh cells ----

gg1 <- rankNet(cellchat_All, mode = "comparison", measure = "weight", title = 'CD4+ macrophage signaling to Tfh cells', stacked = T, do.stat = TRUE, color.use = c('blue', 'red'), targets.use = c('T_CD4_follicular_5'), sources.use = c('Macrophage_CD4neg_28'))
gg2 <- rankNet(cellchat_All, mode = "comparison", measure = "weight", title = 'CD4- macrophage signaling to Tfh cells', stacked = T, do.stat = TRUE, color.use = c('blue', 'red'), targets.use = c('T_CD4_follicular_5'), sources.use = c('Macrophage_CD4pos_25'))
gg3 <- rankNet(cellchat_All, mode = "comparison", measure = "weight", title = 'cDC signaling to Tfh cells', stacked = T, do.stat = TRUE, color.use = c('blue', 'red'), targets.use = c('T_CD4_follicular_5'), sources.use = c('cDC_18'))
gg1+gg2+gg3 

## B cells --> Tfh cells ----
gg1 <- rankNet(cellchat_All, mode = "comparison", measure = "weight", title = 'Resting B signaling to Tfh cells', stacked = T, do.stat = TRUE, color.use = c('blue', 'red'), targets.use = c('T_CD4_follicular_5'), sources.use = c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19'))
gg2 <- rankNet(cellchat_All, mode = "comparison", measure = "weight", title = 'Non-GC cycling B signaling to Tfh cells', stacked = T, do.stat = TRUE, color.use = c('blue', 'red'), targets.use = c('T_CD4_follicular_5'), sources.use = c('B_nonGC_cycling_39'))
gg3 <- rankNet(cellchat_All, mode = "comparison", measure = "weight", title = 'GC LZ B signaling to Tfh cells', stacked = T, do.stat = TRUE, color.use = c('blue', 'red'), targets.use = c('T_CD4_follicular_5'), sources.use = c('B_GC_LZ_2', 'B_GC_LZ_4'))
gg4 <- rankNet(cellchat_All, mode = "comparison", measure = "weight", title = 'GC DZ B macrophage signaling to Tfh cells', stacked = T, do.stat = TRUE, color.use = c('blue', 'red'), targets.use = c('T_CD4_follicular_5'), sources.use = c('B_GC_DZ_1', 'B_GC_DZ_9'))
gg1+gg2+gg3+gg4

### View session information ----
report(sessionInfo())
