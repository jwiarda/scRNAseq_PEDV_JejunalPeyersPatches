library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(scales)
library(report)
library(future)

# Load Seurat object ----
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples.h5Seurat')

# Cluster cells ----

DefaultAssay(seu) <- 'SCT'
seu <- FindNeighbors(seu, reduction = "integrated.cca", dims = 1:30)
plan("multisession", workers = 4)
options(future.globals.maxSize = 30000 * 1024^2, future.seed=TRUE)
seu <- FindClusters(seu, resolution = c(1.5), cluster.name = "cca_clusters_res1.5", graph.name = 'SCT_snn')
DimPlot(seu, group.by = 'cca_clusters_res1.5', reduction = 'umap.cca', label = TRUE, shuffle = TRUE)

# Cluster canonical gene query ----
DefaultAssay(seu) <- 'SCT'
Idents(seu) <- seu$cca_clusters_res1.5
DotPlot(seu, 
        features = c('PTPRC',
                     'CD79A', 'CD79B', 'CD19', 'MS4A1', 'PAX5', # B cells
                     'JCHAIN', 'XBP1', 'PRDM1', 'IRF4', # ASCs
                     'CD69', 'CD83', 'SLA-DQB1', 'SLA-DRA', # B activation
                     'GPR183', 'CCR7', 'KLF2', 'SELL', 'FCER2', 'CD40', # resting B
                     'AICDA', 'CD86', 'BCL6', # follicular B
                     'PCLAF', 'BIRC5', 'TOP2A', 'STMN1'), # cycling B
        cols = c('gold', 'darkmagenta')) + RotatedAxis() 

DotPlot(seu,
        features = c('CD3E', 'CD3G', 'CD247', # T cell
                     'CD4', 'CD8B', 'TRDC', 'CD2', 'CD8A', # T cell subsets
                     'PCLAF', 'BIRC5', 'TOP2A', 'STMN1', # cycling
                     'CCL5', 'ITGAE', # effector/resident
                     'GZMA-16903', 'GZMB', 'GNLY', # cytotoxic
                     'CTSW', 'XCL1', 'SLA-DRA', 'SLA-DQB1', 'CCR9', 'KLRK1', # activation
                     'FCER1G', 'KLRG1', 'ITGB1', 'ITGB7', 'SELL', # SELLhi gd
                     'ID3', 'RHEX', 'BLK', 'SAMSN1', 'IL26', # CD2- gd
                     'CCR7', 'S1PR1', 'LEF1', 'KLF2', # naive ab
                     'BCL6', 'ICOS', 'CTLA4', 'CD40LG', 'IL10', # Non-naive + follicular CD4
                     'CD52', 'IFITM3', 'GPR183', # non-naive CD4
                     'PDCD1', 'CXCR4', 'CD69', # follicular CD4
                     'LTB', 'ID2', 'KIT', 'IL7R', 'IL22', 'KLRB1', 'RORC', 'CXCL8'), # group 3 ILC
        cols = c('gold', 'darkmagenta')) + RotatedAxis()

DotPlot(seu, 
        features = c('AIF1', 'CST3', 'SIRPA',  #myeloid
                     'FLT3', 'SLA-DRA', 'SLA-DQB1', # DC
                     'CD68', 'CXCL2', 'CD4', 'C1QA', 'C1QB', 'C1QC', # macrophage
                     'TCF4', 'XBP1', 'CLEC12A', 'CD93', 'IRF8', 'CD8B', #pDC
                     'ICAM1', 'CSF2RB', 'MS4A2', 'FCER1A'), # mast
        cols = c('gold', 'darkmagenta')) + RotatedAxis()

DotPlot(seu, 
        features = c('PTPRC',
                     'EEF1B2', 'PIGR', 'LYZ', # crypt
                     'LGR5', 'OLFM4', # stem cell
                     'EPCAM', 'FABP2', 'FABP1', 'CLCA4', 'SLC5A1', 'SI', 'ACE2', # enterocyte
                     'GUCA2A', 'GUCA2B', 'BEST4', 'CFTR', 'NOTCH2', 'OTOP2', # BEST4+ enterocyte
                     'TFF3', 'REG4', 'CLCA1', 'SPINK4', 'MUC2', 'CXCL8', 'ITLN2', #goblet
                     'PYY', 'GAST', 'SST', 'CCK','TTR', 'NTS', # EE
                     'NEUROD1', 'CHGA', 'CHGB', 'KRT7', 'SCT', 'PENK'), # EE
        cols = c('gold', 'darkmagenta')) + RotatedAxis()

DotPlot(seu, 
        features = c('PECAM1', 'CDH5', # endothelial
                     'ECM1', 'COL1A1', #fibroblast
                     'MYH11', 'ACTG2'), # muscle
        cols = c('gold', 'darkmagenta')) + RotatedAxis()

# Assign cell annotations ----
seu$celltype <- seu$cca_clusters_res1.5
Idents(seu) <- seu$celltype
types <- c('B_resting' ,'B_GC_DZ', 'B_GC_LZ', 'Enterocyte_mature', 'B_GC_LZ', 'T_CD4_follicular',
           'T_CD4_nonnaive', 'B_resting', 'T_CD8ab', 'B_GC_DZ', 'Enterocyte_intermediate',
           'ILC_group1_ITGAEpos', 'T_ab_resting', 'ILC_group3', 'B_resting', 'B_resting', 
           'T_gd_CD2pos', 'Enterocyte_early', 'cDC', 'B_resting', 'T_gd_CD2pos',
           'ASC', 'ISC_TA', 'ILC_group1_ITGAEpos', 'T_CD8ab', 'Macrophage_CD4pos',
           'T_gd_CD2pos_SELLpos', 'ASC', 'Macrophage_CD4neg', 'T_CD4_cycling', 'Goblet', 
           'Fibroblast', 'B_lowFeature', 'Enterocyte_BEST4pos', 'Goblet', 'Endothelial', 
           'T_gd_CD2neg', 'ILC_group1_ITGAEneg', 'Enterocyte_mature', 'B_nonGC_cycling', 'ASC',
           'ASC', 'Mast') # Rename clusters based on phenotype IDs
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$celltype <- Idents(seu)
Idents(seu) <- seu$celltype
levels(seu) <- c('B_resting', 'B_nonGC_cycling', 'B_GC_LZ', 'B_GC_DZ', 'ASC', 'B_lowFeature',
                 'T_ab_resting', 'T_CD4_nonnaive', 'T_CD4_follicular', 'T_CD4_cycling',
                 'T_CD8ab', 'T_gd_CD2pos', 'T_gd_CD2pos_SELLpos', 'T_gd_CD2neg',
                 'ILC_group1_ITGAEpos', 'ILC_group1_ITGAEneg', 'ILC_group3', 
                 'Macrophage_CD4pos', 'Macrophage_CD4neg', 'cDC', 'Mast',
                 'ISC_TA', 'Enterocyte_early', 'Enterocyte_intermediate', 'Enterocyte_mature', 'Enterocyte_BEST4pos', 'Goblet',
                 'Endothelial', 'Fibroblast')
seu$celltype <- Idents(seu)
DimPlot(seu, group.by = 'celltype', reduction = 'umap.cca', label = TRUE, shuffle = TRUE, 
        cols = c('gold3', 'deeppink4', 'chartreuse4', 'cyan4', 'sandybrown',
                 'cornflowerblue', 'navy', 'lightpink', 'salmon', 
                 'deepskyblue2', 'tan4', 'mediumpurple1', 
                 'darkgreen', 'gray50', 'darkmagenta', 'red', 'hotpink', 'khaki', 
                 'limegreen', 'cadetblue3', 'firebrick', 'deepskyblue4', 'darkseagreen', 'burlywood3', 
                 'goldenrod3', 'blue', 'deeppink', 'lightskyblue3', 'mistyrose3', 'purple'))

# Assign cell lineages ----
seu$celllineage <- seu$celltype
Idents(seu) <- seu$celllineage
types <- c(rep('B/ASC', 6),
           rep('T/ILC', 11),
           rep('Myeloid', 4),
           rep('Epithelial', 6),
           rep('Stromal', 2)) # Rename clusters based on phenotype IDs
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$celllineage <- Idents(seu)
Idents(seu) <- seu$celllineage
levels(seu) <- c('B/ASC', 'T/ILC', 'Myeloid', 'Epithelial', 'Stromal')
seu$celllineage <- Idents(seu)
DimPlot(seu, group.by = 'celllineage', reduction = 'umap.cca', label = TRUE, shuffle = TRUE, 
        cols = c('mediumorchid', 'orange', 'blue', 'chartreuse3', 'grey'))

# Plot sample IDs
DimPlot(seu, group.by = 'AnimalID', reduction = 'umap.cca', label = TRUE, shuffle = TRUE, 
        cols = c('lightblue3', 'steelblue3', 'blue', 'rosybrown2', 'salmon', 'peru', 'red', 'red4'))

# Cell type canonical gene query ----

## Assign new cluster order ----
seu$clusters_celltype <- seu$cca_clusters_res1.5
Idents(seu) <- seu$clusters_celltype
levels(seu) <- c('0', '7', '14', '15', '19', # B_resting
                 '39', # B_nonGC_cycling
                 '2', '4', # B_GC_LZ
                 '1', '9', # B_GC_DZ
                 '21', '27', '40', '41', # ASC
                 '32', # B_lowFeature
                 '12', # T_ab_resting
                 '6', # T_CD4_nonnaive
                 '5', # T_CD4_follicular
                 '29', # T_CD4_cycling
                 '8', '24', # T_CD8ab
                 '16', '20', # T_gd_CD2pos
                 '26', # T_gd_CD2pos_SELLpos
                 '36', # T_gd_CD2neg
                 '11', '23', # ILC_group1_ITGAEpos
                 '37', # ILC_group1_ITGAEneg
                 '13', # ILC_group3
                 '25', # Macrophage_CD4pos
                 '28', # Macrophage_CD4neg
                 '18', # cDC
                 '42', # Mast
                 '22', # ISC_TA
                 '17', # Enterocyte_early
                 '10', # Enterocyte_intermediate
                 '3', '38', # Enterocyte_mature
                 '33', # Enterocyte_BEST4pos
                 '30', '34', # Goblet
                 '35', # Endothelial
                 '31') # Fibroblast
seu$clusters_celltype <- Idents(seu)
Idents(seu) <- seu$clusters_celltype
types <- c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19', # B_resting
           'B_nonGC_cycling_39', # B_nonGC_cycling
           'B_GC_LZ_2', 'B_GC_LZ_4', # B_GC_LZ
           'B_GC_DZ_1', 'B_GC_DZ_9', # B_GC_DZ
           'ASC_21', 'ASC_27', 'ASC_40', 'ASC_41', # ASC
           'B_lowFeature_32', # B_lowFeature
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
           'Fibroblast_31')
names(types) <- levels(seu) # assign GutCellTypes to cluster numbers
seu <- RenameIdents(seu, types) # change dataset identity to cell types in new Seurat object
seu$clusters_celltype <- Idents(seu)
Idents(seu) <- seu$clusters_celltype
levels(seu) <- c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19', # B_resting
                 'B_nonGC_cycling_39', # B_nonGC_cycling
                 'B_GC_LZ_2', 'B_GC_LZ_4', # B_GC_LZ
                 'B_GC_DZ_1', 'B_GC_DZ_9', # B_GC_DZ
                 'ASC_21', 'ASC_27', 'ASC_40', 'ASC_41', # ASC
                 'B_lowFeature_32', # B_lowFeature
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
                 'Fibroblast_31')
seu$clusters_celltype <- Idents(seu)

DefaultAssay(seu) <- 'SCT'
Idents(seu) <- seu$clusters_celltype
DotPlot(seu, 
        features = unique(c('PTPRC',
                            'CD79A', 'CD79B', 'CD19', 'MS4A1', 'PAX5', # B cells
                            'JCHAIN', 'XBP1', 'PRDM1', 'IRF4', # ASCs
                            'CD69', 'CD83', 'SLA-DQB1', 'SLA-DRA', # B activation
                            'GPR183', 'CCR7', 'KLF2', 'SELL', 'FCER2', 'CD40', # resting B
                            'AICDA', 'CD86', 'BCL6', # follicular B
                            'PCLAF', 'BIRC5', 'TOP2A', 'STMN1', # cycling B
                            
                            'CD3E', 'CD3G', 'CD247', # T cell
                            'CD4', 'CD8B', 'TRDC', 'CD2', 'CD8A', # T cell subsets
                            'PCLAF', 'BIRC5', 'TOP2A', 'STMN1', # cycling
                            'CCL5', 'ITGAE', # effector/resident
                            'GZMA-16903', 'GZMB', 'GNLY', # cytotoxic
                            'CTSW', 'XCL1', 'SLA-DRA', 'SLA-DQB1', 'CCR9', 'KLRK1', # activation
                            'FCER1G', 'KLRG1', 'ITGB1', 'ITGB7', 'SELL', # SELLhi gd
                            'ID3', 'RHEX', 'BLK', 'SAMSN1', 'IL26', # CD2- gd
                            'CCR7', 'S1PR1', 'LEF1', 'KLF2', # naive ab
                            'ICOS', 'CTLA4', 'CD40LG', 'IL10', # Non-naive + follicular CD4
                            'CD52', 'IFITM3', 'GPR183', # non-naive CD4
                            'PDCD1', 'CXCR4', 'CD69', # follicular CD4
                            'LTB', 'ID2', 'KIT', 'IL7R', 'IL22', 'KLRB1', 'RORC', 'CXCL8', # group 3 ILC
                            
                            'FLT3', 'SLA-DRA', 'SLA-DQB1', # DC
                            'SIRPA', 'CD68', 'CXCL2', 'C1QA', 'C1QB', 'C1QC', # macrophage
                            'ICAM1', 'CSF2RB', 'MS4A2', 'FCER1A', # mast
                            
                            'RPL5', 'RPS6', 'EEF1B2', 'OLFM4', 'PIGR', 'LYZ', # crypts
                            'FABP2', 'FABP1', 'CLCA4', 'SLC5A1', 'SI', 'ACE2',  # enterocyte
                            'GUCA2A', 'GUCA2B', 'BEST4', 'CFTR', 'NOTCH2', 'OTOP2', # best4 enterocyte
                            'TFF3', 'REG4', 'CLCA1', 'SPINK4', 'MUC2', 'CXCL8', # goblet
                            'PYY', 'GAST', 'SST', 'CCK', 'TTR', 'NTS', # NEUROD1lo EE
                            'NEUROD1', 'CHGA', 'CHGB', 'KRT7', 'SCT', 'PENK', # NEUROD1hi EE
                            
                            'PECAM1', 'CDH5', # endothelial
                            'ECM1', 'COL1A1', 'COL1A2', # fibroblast
                            'TAGLN', 'MYH11', 'ACTG2')), # muscle
        cols = c('gold', 'darkmagenta')) + RotatedAxis() 

# Save file ----
SaveH5Seurat(seu, '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples_Annotated.h5Seurat', overwrite = TRUE)

### View session information ----
report(sessionInfo())