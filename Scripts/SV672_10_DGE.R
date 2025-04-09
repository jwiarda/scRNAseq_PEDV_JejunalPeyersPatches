library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(report)
library(future)
library(presto) # do DGE analysis way fast!
library(stringr)
library(scales)

# Load Seurat object ----
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples_Annotated.h5Seurat')
DefaultAssay(seu) <- 'SCT'

# Cluster DGE ----
## Global ----
de <- wilcoxauc(seu, group_by = 'clusters_celltype', seurat_assay = 'SCT', assay = 'data')
#de <- subset(de, padj < 0.05 & abs(logFC) > 0.25)
write_xlsx(de, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_Global_Clusters_unflitered.xlsx')

## B lineage ----
#prop.table(table(seu$celllineage, seu$clusters_celltype))
de <- wilcoxauc(seu, group_by = 'clusters_celltype', seurat_assay = 'SCT', assay = 'data', 
                groups_use = c('0', '1', '2', '4', '7', '9', '14', '15', '19', '21', '27', '32', '39', '40'))
#de <- subset(de, padj < 0.05 & abs(logFC) > 0.25)
write_xlsx(de, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_Global_Clusters_B_unflitered.xlsx')

## T/ILC lineage ----
#prop.table(table(seu$celllineage, seu$clusters_celltype))
de <- wilcoxauc(seu, group_by = 'clusters_celltype', seurat_assay = 'SCT', assay = 'data', 
                groups_use = c('5', '6', '8', '11', '12', '13', '16', '20', '23', '24', '26', '29', '36', '37'))
#de <- subset(de, padj < 0.05 & abs(logFC) > 0.25)
write_xlsx(de, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_Global_Clusters_TILC_unflitered.xlsx')

## Myeloid lineage ----
#prop.table(table(seu$celllineage, seu$clusters_celltype))
de <- wilcoxauc(seu, group_by = 'clusters_celltype', seurat_assay = 'SCT', assay = 'data', 
                groups_use = c('18', '25', '28', '42'))
#de <- subset(de, padj < 0.05 & abs(logFC) > 0.25)
write_xlsx(de, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_Global_Clusters_Myeloid_unflitered.xlsx')

## Epithelial lineage ----
#prop.table(table(seu$celllineage, seu$clusters_celltype))
de <- wilcoxauc(seu, group_by = 'clusters_celltype', seurat_assay = 'SCT', assay = 'data', 
                groups_use = c('3', '10', '17', '22', '30', '33', '34', '38'))
#de <- subset(de, padj < 0.05 & abs(logFC) > 0.25)
write_xlsx(de, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_Global_Clusters_Epithelial_unflitered.xlsx')

## Stromal lineage ----
#prop.table(table(seu$celllineage, seu$clusters_celltype))
de <- wilcoxauc(seu, group_by = 'clusters_celltype', seurat_assay = 'SCT', assay = 'data', 
                groups_use = c('31', '35'))
#de <- subset(de, padj < 0.05 & abs(logFC) > 0.25)
write_xlsx(de, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_Global_Clusters_Stromal_unflitered.xlsx')

## Pairwise ----
Idents(seu) <- seu$clusters_celltype
clusters <- unique(Idents(seu)) # identify all of our cluster IDs
pairwise <- combn(clusters, 2) # create all pairwise combinations of cluster IDs
p1 <- pairwise[1,] 
p2 <- pairwise[2,] 
comps1 <- data.frame(p1, p2)
colnames(comps1) <- c('pop1', 'pop2')
comps2 <- data.frame(p2, p1)
colnames(comps2) <- c('pop1', 'pop2')
comps <- rbind(comps1, comps2)
dim(comps) # see how many comparisons are to be made, then break them up so we don't have any session time out issues

results <- list()
for(i in 1:nrow(comps)) {
  markers <- wilcoxauc(seu, group_by = 'clusters_celltype', seurat_assay = 'SCT', assay = 'data', groups_use = c(comps[i,1], comps[i,2]))
  markers$pop1 <- paste(comps[i,1])
  markers$pop2 <- paste(comps[i,2])
  results[[i]] <- markers
} 
pwAll <- do.call(rbind, results)
pwAll <- subset(pwAll, logFC > 0) # only take positive FC values because we also perform reciprical comparisons alreday
write.table(pwAll, file='/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_Global_ClustersPairwise_unflitered.tsv', quote=FALSE, sep='\t', col.names = NA)
#write_xlsx(pwAll, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_Global_ClustersPairwise_unflitered.xlsx') # too large to make Excel file
rm(markers, results, pwAll, comps, comps1, comps2, pairwise, clusters)

## Treatment-based DGE within each cluster ----
seu$Infection <- paste(seu$Treatment, seu$LocalInfection_IHC, sep = '_')
Idents(seu) <- seu$Infection
seu <- subset(seu, idents = c('Mock_Uninfected', 'PEDV_PEDV'))
#table(seu$AnimalID) # ensure animal 924 gone
seu$combo <- paste(seu$clusters_celltype, seu$Treatment, sep = '_')
comps <- data.frame(
  (paste(unique(seu$clusters_celltype), 'PEDV', sep = '_')), 
  (paste(unique(seu$clusters_celltype), 'Mock', sep = '_')))

results <- list()
for(i in 1:nrow(comps)) {
  markers <- wilcoxauc(seu, group_by = 'combo', seurat_assay = 'SCT', assay = 'data', groups_use = c(comps[i,1], comps[i,2]))
  markers$pop1 <- paste(comps[i,1])
  markers$pop2 <- paste(comps[i,2])
  results[[i]] <- markers
} 
pwAll <- do.call(rbind, results)
pwAll <- subset(pwAll, logFC > 0) # only take positive FC values because we also perform reciprical comparisons alreday
write.table(pwAll, file='/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_TreatmentWithinClusters_unflitered.tsv', quote=FALSE, sep='\t', col.names = NA)
write_xlsx(pwAll, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_TreatmentWithinClusters_unflitered.xlsx')
#seu$combo <- NULL

### Heatmap of treatment-based DGE within each cluster ----

de <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_TreatmentWithinClusters_unflitered.tsv')
de <- subset(de, abs(logFC) > 0.25 & padj < 0.05)
de <- subset(de, pct_in > 0.1 | pct_out > 0.1)
de$group <- sub("_[^_]*$", "", de$group)
de <- subset(de, subset = !(group == 'B_lowFeature_32'))
de <- data.frame(table(de$group))
de$Var1 <- factor(de$Var1, levels = rev(c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19', # B_resting
                                          'B_nonGC_cycling_39', # B_nonGC_cycling
                                          'B_GC_LZ_2', 'B_GC_LZ_4', # B_GC_LZ
                                          'B_GC_DZ_1', 'B_GC_DZ_9', # B_GC_DZ
                                          'ASC_21', 'ASC_27', 'ASC_40', 'ASC_41', # ASC
                                          #'B_lowFeature_32', # B_lowFeature
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

ggplot(de, aes(y = Var1, x = 1, fill = Freq)) +
  geom_tile(color = 'black')+
  geom_text(aes(label = Freq)) +
  scale_fill_gradientn('# of DEGs', colours = c('yellow', 'orange', 'red', 'red3'),
                       limits = c(0, 1500), oob=squish)+ 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold",
                                   angle=0, vjust=.5, hjust=1),
        axis.title.y = element_text(size = 20)) +
  labs(y = "Cluster") +
  ggtitle('Mock vs PEDV DEGs')

### View session information ----
report(sessionInfo())
