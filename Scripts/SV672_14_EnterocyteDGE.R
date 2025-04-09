library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(report)
library(presto) # do DGE analysis way fast!
library(eulerr)

# Load Seurat object ----
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/trajectorySubset_PolarizedEnterocytes.h5seurat')
DefaultAssay(seu) <- 'SCT'

# Perform DGE across polarized states ----
# Define polarized states
seu$pseudotime1[is.na(seu$pseudotime1)] <- -5
seu$pseudotime2[is.na(seu$pseudotime2)] <- -5
seu$EnterocytePolarization <- 'Both'
seu$EnterocytePolarization[seu$pseudotime1 > 0 & seu$pseudotime2 < 0] <- 'Traj1'
seu$EnterocytePolarization[seu$pseudotime2 > 0 & seu$pseudotime1 < 0] <- 'Traj2'

# DGE
de <- wilcoxauc(seu, group_by = 'EnterocytePolarization', seurat_assay = 'SCT', assay = 'data', groups_use = c('Traj1', 'Traj2'))
write.table(de, file='/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_Global_PolarEnterocyteTrajectoryStatesPairwise_unflitered.tsv', quote=FALSE, sep='\t', col.names = NA)
de <- subset(de, logFC > 0.25 & padj < 0.05)
de <- subset(de, pct_in > 0.1 | pct_out > 0.1)

## Plot gene & mitochondrial expression ----

# % mito
annotKey <- read_excel('/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx') # read in our key for modified gene annotations
annot <- "/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS_PEDV_USA_NC_2013_49469.gtf" # specify file path to Sus scrofa 11.1 version 97 annotation file with modified gene names
mitoGenes <-  system2("grep", args = c('^MT', annot, "| grep -o 'ENSSSCG[0-9]*' | uniq"), stdout = TRUE) # extract mitochondrial gene Ensembl IDs from annotation file
mitoGenes <- annotKey[annotKey$ENSID %in% mitoGenes,] # identify rows for mitochondrial genes in the annotation key
mitoGenes <- mitoGenes$FinalList
length(mitoGenes) # make sure the length is 37, the same number of mitochondrial genes found in pigs as in humans
rm(annot, annotKey) # clear space

# Now we need to calculate the percentage of mitochondrial read counts in each cell:
counts <- GetAssayData(object = seu, layer = "counts")
mitoCounts <- counts[rownames(counts) %in% mitoGenes,] # identify rows for mitochondrial genes in the annotation key
pctMito <- ((colSums(mitoCounts))/(colSums(counts)))*100
seu <- AddMetaData(seu, pctMito, col.name = "percent_mito")
rm(counts, mitoCounts,pctMito,mitoGenes)

# top DE genes for Traj1
Idents(seu) <- seu$EnterocytePolarization
sub <- subset(seu, idents = c('Traj1', 'Traj2'))
sub <- ScaleData(sub)
topgenes <- de %>% group_by(group) %>% arrange(group, desc(logFC)) %>% top_n(10, logFC) 
topgenes <- subset(topgenes, group == 'Traj1')

# Plot
sub$EnterocyteState <- sub$EnterocytePolarization
sub$EnterocyteState[sub$EnterocytePolarization == 'Traj1'] <- 'HomeostaticEnterocyte'
sub$EnterocyteState[sub$EnterocytePolarization == 'Traj2'] <- 'StressedEnterocyte'
Idents(sub) <- sub$EnterocyteState
VlnPlot(sub, features = c(topgenes$feature, 'percent_mito'), cols = c('limegreen', 'purple1'), ncol = 5)

# Perform DGE across infection states ----
# Define infection states
sub$Infection <- NA
sub$Infection[sub$Treatment == 'PEDV' & sub[["SCT"]]@data['PEDV',] > 0] <- 'Infected'
sub$Infection[sub$Treatment == 'PEDV' & sub[["SCT"]]@data['PEDV',] == 0] <- 'Bystander'
sub$Infection[sub$Treatment == 'Mock'] <- 'Mock'

sub$EnterocyteState <- sub$EnterocytePolarization
sub$EnterocyteState[sub$EnterocytePolarization == 'Traj1'] <- 'HomeostaticEnterocyte'
sub$EnterocyteState[sub$EnterocytePolarization == 'Traj2'] <- 'StressedEnterocyte'

sub$EnterocyteInfection <- paste(sub$EnterocyteState, sub$Infection, sep = '_')

# DGE

clusters <- c('HomeostaticEnterocyte_Mock', 'HomeostaticEnterocyte_Bystander', 'HomeostaticEnterocyte_Infected')
comps <- data.frame(t(combn(clusters, 2))) # create all pairwise combinations of cluster IDs
colnames(comps) <- c('pop1', 'pop2')

results <- list()
for(i in 1:nrow(comps)) {
  markers <- wilcoxauc(sub, group_by = 'EnterocyteInfection', seurat_assay = 'SCT', assay = 'data', groups_use = c(comps[i,1], comps[i,2]))
  markers$pop1 <- paste(comps[i,1])
  markers$pop2 <- paste(comps[i,2])
  results[[i]] <- markers
} 
hs <- do.call(rbind, results)

clusters <- c('StressedEnterocyte_Mock', 'StressedEnterocyte_Bystander', 'StressedEnterocyte_Infected')
comps <- data.frame(t(combn(clusters, 2))) # create all pairwise combinations of cluster IDs
colnames(comps) <- c('pop1', 'pop2')

results <- list()
for(i in 1:nrow(comps)) {
  markers <- wilcoxauc(sub, group_by = 'EnterocyteInfection', seurat_assay = 'SCT', assay = 'data', groups_use = c(comps[i,1], comps[i,2]))
  markers$pop1 <- paste(comps[i,1])
  markers$pop2 <- paste(comps[i,2])
  results[[i]] <- markers
} 
st <- do.call(rbind, results)

pwAll <- rbind(hs, st)
write.table(pwAll, file='/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_ClustersPairwise_PEDVEnterocyteInfectionStates_unflitered.tsv', quote=FALSE, sep='\t', col.names = NA)
pwAll <- subset(pwAll, logFC > 0.25 & padj < 0.05)
pwAll <- subset(pwAll, pct_in > 0.1 | pct_out > 0.1)

# Venn Diagram ----

Hs_BystanderVsMock <- subset(pwAll, group == 'HomeostaticEnterocyte_Bystander' & (pop1 == 'HomeostaticEnterocyte_Mock' | pop2 == 'HomeostaticEnterocyte_Mock'))
#Hs_BystanderVsInfected <- subset(pwAll, group == 'HomeostaticEnterocyte_Bystander' & (pop1 == 'HomeostaticEnterocyte_Infected' | pop2 == 'HomeostaticEnterocyte_Infected'))
Hs_InfectedVsMock <- subset(pwAll, group == 'HomeostaticEnterocyte_Infected' & (pop1 == 'HomeostaticEnterocyte_Mock' | pop2 == 'HomeostaticEnterocyte_Mock'))
#Hs_InfectedVsBystander <- subset(pwAll, group == 'HomeostaticEnterocyte_Infected' & (pop1 == 'HomeostaticEnterocyte_Bystander' | pop2 == 'HomeostaticEnterocyte_Bystander'))
#Hs_MockVsBystander <- subset(pwAll, group == 'HomeostaticEnterocyte_Mock' & (pop1 == 'HomeostaticEnterocyte_Bystander' | pop2 == 'HomeostaticEnterocyte_Bystander'))
#Hs_MockVsInfected <- subset(pwAll, group == 'HomeostaticEnterocyte_Mock' & (pop1 == 'HomeostaticEnterocyte_Infected' | pop2 == 'HomeostaticEnterocyte_Infected'))

St_BystanderVsMock <- subset(pwAll, group == 'StressedEnterocyte_Bystander' & (pop1 == 'StressedEnterocyte_Mock' | pop2 == 'StressedEnterocyte_Mock'))
#St_BystanderVsInfected <- subset(pwAll, group == 'StressedEnterocyte_Bystander' & (pop1 == 'StressedEnterocyte_Infected' | pop2 == 'StressedEnterocyte_Infected'))
St_InfectedVsMock <- subset(pwAll, group == 'StressedEnterocyte_Infected' & (pop1 == 'StressedEnterocyte_Mock' | pop2 == 'StressedEnterocyte_Mock'))
#St_InfectedVsBystander <- subset(pwAll, group == 'StressedEnterocyte_Infected' & (pop1 == 'StressedEnterocyte_Bystander' | pop2 == 'StressedEnterocyte_Bystander'))
#St_MockVsBystander <- subset(pwAll, group == 'StressedEnterocyte_Mock' & (pop1 == 'StressedEnterocyte_Bystander' | pop2 == 'StressedEnterocyte_Bystander'))
#St_MockVsInfected <- subset(pwAll, group == 'StressedEnterocyte_Mock' & (pop1 == 'StressedEnterocyte_Infected' | pop2 == 'StressedEnterocyte_Infected'))

l <- tibble::lst(unique(Hs_BystanderVsMock$feature), unique(Hs_InfectedVsMock$feature),
                 unique(St_BystanderVsMock$feature), unique(St_InfectedVsMock$feature))
dat <- data.frame(lapply(l, `length<-`, max(lengths(l))))
dat <- as.data.frame.matrix(table(stack(dat)))
colnames(dat) <- c('Homeostatic Bystander', 'Homeostatic Infected', 'Stressed Bystander', 'Stressed Infected')
#plot(eulerr::euler(dat), quantities = TRUE, alpha = 0.3, fill = c('dodgerblue', 'forestgreen', 'grey70', 'darkmagenta'), alpha = 0.3) 
plot(eulerr::venn(dat), quantities = TRUE, alpha = 0.3, fill = c('dodgerblue', 'forestgreen', 'grey70', 'darkmagenta'), alpha = 0.3) 

### View session information ----
report(sessionInfo())