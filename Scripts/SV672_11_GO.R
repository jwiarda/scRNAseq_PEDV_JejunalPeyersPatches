library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(report)
library(future)
library(stringr)
library(scales)
library(tidyverse)
library(biomaRt)  # used to map GOterms to Ensembl IDs
library(viridis)

# Load Seurat object ----
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples_Annotated.h5Seurat')
DefaultAssay(seu) <- 'SCT'

## Read in & filter DE data ----
AllCellsDGE <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_TreatmentWithinClusters_unflitered.tsv')
AllCellsDGE <- subset(AllCellsDGE, abs(logFC) > 0.25 & padj < 0.05)
AllCellsDGE <- subset(AllCellsDGE, pct_in > 0.1)

## GO analysis ----
dir.create('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/')

GO_background_generator <- 
  function(GO_term_universe, seurat_obj, output_path){
    # gene list is enriched genes, 
    # GO_term_universe is already prepared, needs to be filtered
    # seurat obj will be used to filter GO_term_universe
    # browser()
    # need to make sure all genes in seurat object are detected
    counts <- seurat_obj@assays$SCT@counts
    if (!all(rowSums(counts) > 0)){
      print('some genes not detected, filtering genes with 0 counts')
      counts <- counts[rowSums(counts) > 0,]
    }
    
    all_GO <- read_tsv('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_gene_to_GO.tsv')
    new_universe <-  all_GO %>% 
      filter(FinalList %in% rownames(counts)) %>% 
      write_tsv(output_path)
    print(paste('number of genes in all universe',nrow(all_GO)))
    print(paste('number of genes in new universe',nrow(new_universe)))
    
  }

# wrapper function that makes using topGO easier (I think)

# you can also access this function by loading my R package 'funfuns'
# that way you don't need to read in this function every time
# remotes::install_github('jtrachsel/funfuns')
# library(funfuns)

topGO_wrapper <- function(myInterestingGenes, #vector
                          mapping_file,       # two column file
                          ont='BP',
                          algor = 'elim',
                          statistic='Fisher',
                          nodeSize=10, 
                          return_GOdata=F){
  
  require(topGO)
  # browser()
  geneID2GO <- readMappings(mapping_file)
  geneNames <- names(geneID2GO)
  
  # Get the list of genes of interest
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  
  #initialize topGOdata object
  GOdata <- new("topGOdata", ontology = ont, allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO,
                nodeSize=nodeSize)
  
  # contstruct a tibble that maps interesting genes to GO terms
  # this can add a decent amount of time to the function call...
  interesting_genes_in_GOs <- 
    genesInTerm(GOdata) %>% 
    enframe(name='GO.ID', 
            value='genes_in_GO') %>% 
    mutate(involved_genes=map(.x=genes_in_GO, ~ .x[.x %in% myInterestingGenes]), 
           involved_genes=map_chr(.x=involved_genes, ~paste(.x, collapse = '_'))) %>% 
    dplyr::select(GO.ID, involved_genes) 
  
  # Run topGO test
  resultTopGO.elim <- runTest(GOdata, algorithm = algor, statistic = statistic )
  allRes <- GenTable(GOdata, pval = resultTopGO.elim,
                     orderBy = "pval",
                     topNodes = length(GOdata@graph@nodes), #include all nodes
                     numChar=1000)
  
  # clean up results and add in extra info
  allRes <- allRes %>%
    mutate(ont=ifelse(ont=='BP', 'Biological Process',
                      ifelse(ont=='MF', 'Molecular Function', "Cellular Component"))) %>%
    mutate(GO_aspect = ont,
           algorithm = algor,
           statistic = statistic) %>%
    dplyr::select(-ont) %>% 
    left_join(interesting_genes_in_GOs)
  
  if (return_GOdata == TRUE){
    return(list(allRes, GOdata))
  } else {
    return(allRes)
  }
  
  
}


#### Get ensemble gene_IDs for pig genome ###
# This section uses the internet to map these IDs, so sometimes it is very slow

# select mart and data set
bm <- useMart("ensembl")
bm <- useDataset("sscrofa_gene_ensembl", mart=bm)

# Get ensembl gene ids and GO terms
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','go_id'))

# examine result
head(EG2GO,15)

# Remove blank entries
EG2GO <- EG2GO[EG2GO$go_id != '',]

### format mapping file for use with topGO wrapper function
EG2GO <- 
  EG2GO %>%
  group_by(ensembl_gene_id) %>% 
  summarise(GO=paste(go_id, sep = ' ', collapse = ',')) %>% 
  transmute(EnsemblID=ensembl_gene_id, 
            GO=GO) 

# scRNA-seq all cells ----

# to map gene 'FinalList' annotation to Ensembl gene_IDs
gene_IDs <- read_excel("/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/UpdatedGeneNameListForSus97GTF_06302021_JEW_SKS.xlsx")
gene_IDs$EnsemblID <- gene_IDs$ENSID # make identical to column name in EG2GO

detected_genes <- rownames(seu)

# all GO terms detected in this dataset

GO_gene_universe <- 
  gene_IDs %>% 
  filter(FinalList %in% detected_genes) %>% 
  left_join(EG2GO) %>% 
  filter(!is.na(GO))

GO_gene_universe %>%
  dplyr::select(FinalList, GO) %>% 
  write_tsv('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_gene_to_GO.tsv')


# Cell type-specific GO terms

AllCells_results <- 
  AllCellsDGE %>% 
  filter(logFC > 0) %>% 
  group_by(group) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(feature))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_gene_to_GO.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(group, filt_results) %>% 
  unnest(cols=filt_results)

AllCells_results$Fold_enrichment <- AllCells_results$Significant/AllCells_results$Expected # add in stat for fold enrichment
AllCells_results <- subset(AllCells_results, Significant > 1) # remove GO terms with only one enriched genedim

AllCells_results %>% write_tsv('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_TreatmentWithinCluster_GOresults.tsv')

# Identify top 10 processes per lineage + treatment ----
# Criteria used was GO term had pval < 0.05 in at least half the populations of a lineage, then took the 10 terms with lowest median p value. For stromal lineage, had to be GO term for all (in this case 2) populations instead of half.

## Plot GO terms in PEDV treatment ----
GO <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_TreatmentWithinCluster_GOresults.tsv')
GO$Treatment <- sub('.*\\_', '', GO$group)
GO$Cluster <- sub("_[^_]+$", "", GO$group)
GO <- subset(GO, Treatment == 'PEDV')
GO$pval <- as.numeric(GO$pval)

vec <- c('ISC_TA_22', # ISC_TA
        'Enterocyte_early_17', # Enterocyte_early
        'Enterocyte_intermediate_10', # Enterocyte_intermediate
        'Enterocyte_mature_3', 'Enterocyte_mature_38', # Enterocyte_mature
        'Enterocyte_BEST4pos_33', # Enterocyte_BEST4pos
        'Goblet_30', 'Goblet_34') # Goblet
subGO <- GO[GO$Cluster %in% vec,] 
subCount <- data.frame(table(subGO$GO.ID))
subCount <- subset(subCount, Freq >= (length(vec)/2))
subGO <- subGO[subGO$GO.ID %in% subCount$Var1,]
ag <- aggregate(pval ~ GO.ID, data = subGO, summary)
ag <- slice_min(ag, ag$pval[, 'Median'], n = 10)
eGO_PEDV <- subGO[subGO$GO.ID %in% ag$GO.ID,]

vec <- c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19', # B_resting
'B_nonGC_cycling_39', # B_nonGC_cycling
'B_GC_LZ_2', 'B_GC_LZ_4', # B_GC_LZ
'B_GC_DZ_1', 'B_GC_DZ_9', # B_GC_DZ
'ASC_21', 'ASC_27', 'ASC_40', 'ASC_41') # ASC
subGO <- GO[GO$Cluster %in% vec,] 
subCount <- data.frame(table(subGO$GO.ID))
subCount <- subset(subCount, Freq >= (length(vec)/2))
subGO <- subGO[subGO$GO.ID %in% subCount$Var1,]
ag <- aggregate(pval ~ GO.ID, data = subGO, summary)
ag <- slice_min(ag, ag$pval[, 'Median'], n = 10)
bGO_PEDV <- subGO[subGO$GO.ID %in% ag$GO.ID,]

vec <- c('T_ab_resting_12', # T_ab_resting
'T_CD4_nonnaive_6', # T_CD4_nonnaive
'T_CD4_follicular_5', # T_CD4_follicular
'T_CD4_cycling_29', # T_CD4_cycling
'T_CD8ab_8', 'T_CD8ab_24', # T_CD8ab
'T_gd_CD2pos_16', 'T_gd_CD2pos_20', # T_gd_CD2pos
'T_gd_CD2pos_SELLpos_26', # T_gd_CD2pos_SELLpos
'T_gd_CD2neg_36', # T_gd_CD2neg
'ILC_group1_ITGAEpos_11', 'ILC_group1_ITGAEpos_23', # ILC_group1_ITGAEpos
'ILC_group1_ITGAEneg_37', # ILC_group1_ITGAEneg
'ILC_group3_13') # ILC_group3
subGO <- GO[GO$Cluster %in% vec,] 
subCount <- data.frame(table(subGO$GO.ID))
subCount <- subset(subCount, Freq >= (length(vec)/2))
subGO <- subGO[subGO$GO.ID %in% subCount$Var1,]
ag <- aggregate(pval ~ GO.ID, data = subGO, summary)
ag <- slice_min(ag, ag$pval[, 'Median'], n = 10)
tGO_PEDV <- subGO[subGO$GO.ID %in% ag$GO.ID,]

vec <- c('Macrophage_CD4pos_25', # Macrophage_CD4pos
'Macrophage_CD4neg_28', # Macrophage_CD4neg
'cDC_18', # cDC
'Mast_42') # Mast
subGO <- GO[GO$Cluster %in% vec,] 
subCount <- data.frame(table(subGO$GO.ID))
subCount <- subset(subCount, Freq >= (length(vec)/2))
subGO <- subGO[subGO$GO.ID %in% subCount$Var1,]
ag <- aggregate(pval ~ GO.ID, data = subGO, summary)
ag <- slice_min(ag, ag$pval[, 'Median'], n = 10)
mGO_PEDV <- subGO[subGO$GO.ID %in% ag$GO.ID,]

vec <- c('Endothelial_35', # Endothelial
'Fibroblast_31')
subGO <- GO[GO$Cluster %in% vec,] 
subCount <- data.frame(table(subGO$GO.ID))
subCount <- subset(subCount, Freq >= (length(vec)))
subGO <- subGO[subGO$GO.ID %in% subCount$Var1,]
ag <- aggregate(pval ~ GO.ID, data = subGO, summary)
ag <- slice_min(ag, ag$pval[, 'Median'], n = 10)
sGO_PEDV <- subGO[subGO$GO.ID %in% ag$GO.ID,]

# Create dot plot ----
## Dot plot just B lineage ----
GO <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_TreatmentWithinCluster_GOresults.tsv')
GO$Treatment <- sub('.*\\_', '', GO$group)
GO$Cluster <- sub("_[^_]+$", "", GO$group)
vec <- c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19', # B_resting
         'B_nonGC_cycling_39', # B_nonGC_cycling
         'B_GC_LZ_2', 'B_GC_LZ_4', # B_GC_LZ
         'B_GC_DZ_1', 'B_GC_DZ_9', # B_GC_DZ
         'ASC_21', 'ASC_27', 'ASC_40', 'ASC_41') # ASC
GO <- GO[GO$Cluster %in% vec,] 
GO$pval <- as.numeric(GO$pval)
GO$Fold_enrichment <- as.numeric(GO$Fold_enrichment)
GO <- GO[GO$GO.ID %in% c(bGO_PEDV$GO.ID),]
table(GO$GO.ID, GO$Treatment)
GO$Process <- paste(GO$GO.ID, GO$Term, sep = ': ')

GO$group <- factor(GO$group, levels = c(paste(c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19', # B_resting
                                                'B_nonGC_cycling_39', # B_nonGC_cycling
                                                'B_GC_LZ_2', 'B_GC_LZ_4', # B_GC_LZ
                                                'B_GC_DZ_1', 'B_GC_DZ_9', # B_GC_DZ
                                                'ASC_21', 'ASC_27', 'ASC_40', 'ASC_41' # ASC
                                                ), 'PEDV', sep = '_'),
                                        paste(c('B_resting_0', 'B_resting_7', 'B_resting_14', 'B_resting_15', 'B_resting_19', # B_resting
                                                'B_nonGC_cycling_39', # B_nonGC_cycling
                                                'B_GC_LZ_2', 'B_GC_LZ_4', # B_GC_LZ
                                                'B_GC_DZ_1', 'B_GC_DZ_9', # B_GC_DZ
                                                'ASC_21', 'ASC_27', 'ASC_40', 'ASC_41' # ASC
                                                ), 'Mock', sep = '_')))
GO$Process <- factor(GO$Process, levels = rev(c("GO:0002474: antigen processing and presentation of peptide antigen via MHC class I",
                                                "GO:0045071: negative regulation of viral genome replication",
                                                "GO:0035455: response to interferon-alpha",
                                                "GO:2000561: regulation of CD4-positive, alpha-beta T cell proliferation",
                                                "GO:0032481: positive regulation of type I interferon production",
                                                "GO:1900151: regulation of nuclear-transcribed mRNA catabolic process, deadenylation-dependent decay",
                                                "GO:0050862: positive regulation of T cell receptor signaling pathway",
                                                "GO:0006955: immune response",
                                                "GO:0045648: positive regulation of erythrocyte differentiation",
                                                "GO:0006511: ubiquitin-dependent protein catabolic process",
                                                "GO:0045577: regulation of B cell differentiation",
                                                "GO:0009968: negative regulation of signal transduction",
                                                "GO:0032733: positive regulation of interleukin-10 production",
                                                "GO:0002699: positive regulation of immune effector process",
                                                "GO:1900225: regulation of NLRP3 inflammasome complex assembly")))

ggplot(GO) +
  geom_point(aes(x = group, 
                 y = Process,
                 size = (1-pval),
                 color = Fold_enrichment)) +
  scale_size(limits=c(.95,1),
             breaks=c(.95, .96,.97,.98,.99,1),
             labels=c(".05",".04",".03",".02",".01","<.01"),
             name = "p-value",
             guide="legend")+
  scale_colour_viridis(option = 'plasma', end = 0.85, begin = 0.2, limits = c(0,50), oob = squish,
                       labels = c('0', '10', '20', '30', '40', '50+'), name = 'Fold enrichment') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  #scale_x_discrete(drop = FALSE)+
  scale_y_discrete(limits=rev)

## Dot plot just T/ILC lineage ----
GO <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_TreatmentWithinCluster_GOresults.tsv')
GO$Treatment <- sub('.*\\_', '', GO$group)
GO$Cluster <- sub("_[^_]+$", "", GO$group)
vec <- c('T_ab_resting_12', # T_ab_resting
         'T_CD4_nonnaive_6', # T_CD4_nonnaive
         'T_CD4_follicular_5', # T_CD4_follicular
         'T_CD4_cycling_29', # T_CD4_cycling
         'T_CD8ab_8', 'T_CD8ab_24', # T_CD8ab
         'T_gd_CD2pos_16', 'T_gd_CD2pos_20', # T_gd_CD2pos
         'T_gd_CD2pos_SELLpos_26', # T_gd_CD2pos_SELLpos
         'T_gd_CD2neg_36', # T_gd_CD2neg
         'ILC_group1_ITGAEpos_11', 'ILC_group1_ITGAEpos_23', # ILC_group1_ITGAEpos
         'ILC_group1_ITGAEneg_37', # ILC_group1_ITGAEneg
         'ILC_group3_13') 
GO <- GO[GO$Cluster %in% vec,] 
GO$pval <- as.numeric(GO$pval)
GO$Fold_enrichment <- as.numeric(GO$Fold_enrichment)
GO <- GO[GO$GO.ID %in% c(tGO_PEDV$GO.ID),]
table(GO$GO.ID, GO$Treatment)
GO$Process <- paste(GO$GO.ID, GO$Term, sep = ': ')

GO$group <- factor(GO$group, levels = c(paste(c('T_ab_resting_12', # T_ab_resting
                                                'T_CD4_nonnaive_6', # T_CD4_nonnaive
                                                'T_CD4_follicular_5', # T_CD4_follicular
                                                'T_CD4_cycling_29', # T_CD4_cycling
                                                'T_CD8ab_8', 'T_CD8ab_24', # T_CD8ab
                                                'T_gd_CD2pos_16', 'T_gd_CD2pos_20', # T_gd_CD2pos
                                                'T_gd_CD2pos_SELLpos_26', # T_gd_CD2pos_SELLpos
                                                'T_gd_CD2neg_36', # T_gd_CD2neg
                                                'ILC_group1_ITGAEpos_11', 'ILC_group1_ITGAEpos_23', # ILC_group1_ITGAEpos
                                                'ILC_group1_ITGAEneg_37', # ILC_group1_ITGAEneg
                                                'ILC_group3_13'
), 'PEDV', sep = '_'),
paste(c('T_ab_resting_12', # T_ab_resting
        'T_CD4_nonnaive_6', # T_CD4_nonnaive
        'T_CD4_follicular_5', # T_CD4_follicular
        'T_CD4_cycling_29', # T_CD4_cycling
        'T_CD8ab_8', 'T_CD8ab_24', # T_CD8ab
        'T_gd_CD2pos_16', 'T_gd_CD2pos_20', # T_gd_CD2pos
        'T_gd_CD2pos_SELLpos_26', # T_gd_CD2pos_SELLpos
        'T_gd_CD2neg_36', # T_gd_CD2neg
        'ILC_group1_ITGAEpos_11', 'ILC_group1_ITGAEpos_23', # ILC_group1_ITGAEpos
        'ILC_group1_ITGAEneg_37', # ILC_group1_ITGAEneg
        'ILC_group3_13'
), 'Mock', sep = '_')))
GO$Process <- factor(GO$Process, levels = rev(c("GO:0045087: innate immune response",
                                            "GO:0140374: antiviral innate immune response",
                                            "GO:0002474: antigen processing and presentation of peptide antigen via MHC class I",
                                            "GO:0001961: positive regulation of cytokine-mediated signaling pathway",
                                            "GO:0032691: negative regulation of interleukin-1 beta production",
                                            "GO:0032760: positive regulation of tumor necrosis factor production",
                                            "GO:0032728: positive regulation of interferon-beta production",
                                            "GO:0060337: type I interferon-mediated signaling pathway",
                                            "GO:0071357: cellular response to type I interferon",
                                            "GO:0035455: response to interferon-alpha",
                                            "GO:0019646: aerobic electron transport chain",
                                            "GO:0042775: mitochondrial ATP synthesis coupled electron transport")))

ggplot(GO) +
  geom_point(aes(x = group, 
                 y = Process,
                 size = (1-pval),
                 color = Fold_enrichment)) +
  scale_size(limits=c(.95,1),
             breaks=c(.95, .96,.97,.98,.99,1),
             labels=c(".05",".04",".03",".02",".01","<.01"),
             name = "p-value",
             guide="legend")+
  scale_colour_viridis(option = 'plasma', end = 0.85, begin = 0.2, limits = c(0,50), oob = squish,
                       labels = c('0', '10', '20', '30', '40', '50+'), name = 'Fold enrichment') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  #scale_x_discrete(drop = FALSE)+
  scale_y_discrete(limits=rev)

## Dot plot just myeloid lineage ----
GO <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_TreatmentWithinCluster_GOresults.tsv')
GO$Treatment <- sub('.*\\_', '', GO$group)
GO$Cluster <- sub("_[^_]+$", "", GO$group)
vec <- c('Macrophage_CD4pos_25', # Macrophage_CD4pos
         'Macrophage_CD4neg_28', # Macrophage_CD4neg
         'cDC_18', # cDC
         'Mast_42') 
GO <- GO[GO$Cluster %in% vec,] 
GO$pval <- as.numeric(GO$pval)
GO$Fold_enrichment <- as.numeric(GO$Fold_enrichment)
GO <- GO[GO$GO.ID %in% c(mGO_PEDV$GO.ID),]
table(GO$GO.ID, GO$Treatment)
GO$Process <- paste(GO$GO.ID, GO$Term, sep = ': ')

GO$group <- factor(GO$group, levels = c(paste(c('Macrophage_CD4pos_25', # Macrophage_CD4pos
                                                'Macrophage_CD4neg_28', # Macrophage_CD4neg
                                                'cDC_18', # cDC
                                                'Mast_42'
), 'PEDV', sep = '_'),
paste(c('Macrophage_CD4pos_25', # Macrophage_CD4pos
        'Macrophage_CD4neg_28', # Macrophage_CD4neg
        'cDC_18', # cDC
        'Mast_42'
), 'Mock', sep = '_')))
GO$Process <- factor(GO$Process, levels = rev(c("GO:0006955: immune response",
                                                "GO:0071353: cellular response to interleukin-4",
                                                "GO:0030593: neutrophil chemotaxis",
                                                "GO:0060337: type I interferon-mediated signaling pathway",
                                                "GO:0009967: positive regulation of signal transduction",
                                                "GO:0007229: integrin-mediated signaling pathway",
                                                "GO:1903077: negative regulation of protein localization to plasma membrane",
                                                "GO:0034314: Arp2/3 complex-mediated actin nucleation",
                                                "GO:2001244: positive regulation of intrinsic apoptotic signaling pathway" ,
                                                "GO:0070936: protein K48-linked ubiquitination")))

ggplot(GO) +
  geom_point(aes(x = group, 
                 y = Process,
                 size = (1-pval),
                 color = Fold_enrichment)) +
  scale_size(limits=c(.95,1),
             breaks=c(.95, .96,.97,.98,.99,1),
             labels=c(".05",".04",".03",".02",".01","<.01"),
             name = "p-value",
             guide="legend")+
  scale_colour_viridis(option = 'plasma', end = 0.85, begin = 0.2, limits = c(0,50), oob = squish,
                       labels = c('0', '10', '20', '30', '40', '50+'), name = 'Fold enrichment') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  #scale_x_discrete(drop = FALSE)+
  scale_y_discrete(limits=rev)

## Dot plot just epithelial lineage ----
GO <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_TreatmentWithinCluster_GOresults.tsv')
GO$Treatment <- sub('.*\\_', '', GO$group)
GO$Cluster <- sub("_[^_]+$", "", GO$group)
vec <- c('ISC_TA_22', # ISC_TA
         'Enterocyte_early_17', # Enterocyte_early
         'Enterocyte_intermediate_10', # Enterocyte_intermediate
         'Enterocyte_mature_3', 'Enterocyte_mature_38', # Enterocyte_mature
         'Enterocyte_BEST4pos_33', # Enterocyte_BEST4pos
         'Goblet_30', 'Goblet_34') 
GO <- GO[GO$Cluster %in% vec,] 
GO$pval <- as.numeric(GO$pval)
GO$Fold_enrichment <- as.numeric(GO$Fold_enrichment)
GO <- GO[GO$GO.ID %in% c(eGO_PEDV$GO.ID),]
table(GO$GO.ID, GO$Treatment)
GO$Process <- paste(GO$GO.ID, GO$Term, sep = ': ')

GO$group <- factor(GO$group, levels = c(paste(c('ISC_TA_22', # ISC_TA
                                                'Enterocyte_early_17', # Enterocyte_early
                                                'Enterocyte_intermediate_10', # Enterocyte_intermediate
                                                'Enterocyte_mature_3', 'Enterocyte_mature_38', # Enterocyte_mature
                                                'Enterocyte_BEST4pos_33', # Enterocyte_BEST4pos
                                                'Goblet_30', 'Goblet_34'
), 'PEDV', sep = '_'),
paste(c('ISC_TA_22', # ISC_TA
        'Enterocyte_early_17', # Enterocyte_early
        'Enterocyte_intermediate_10', # Enterocyte_intermediate
        'Enterocyte_mature_3', 'Enterocyte_mature_38', # Enterocyte_mature
        'Enterocyte_BEST4pos_33', # Enterocyte_BEST4pos
        'Goblet_30', 'Goblet_34'
), 'Mock', sep = '_')))
GO$Process <- factor(GO$Process, levels = rev(c(
                                                "GO:0140374: antiviral innate immune response",
                                                "GO:0006465: signal peptide processing",
                                                "GO:1903078: positive regulation of protein localization to plasma membrane",
                                                "GO:0000028: ribosomal small subunit assembly",
                                                "GO:0006457: protein folding",
                                                "GO:0018279: protein N-linked glycosylation via asparagine",
                                                "GO:0000380: alternative mRNA splicing, via spliceosome",
                                                "GO:0051336: regulation of hydrolase activity",
                                                "GO:1990830: cellular response to leukemia inhibitory factor",
                                                "GO:0007339: binding of sperm to zona pellucida",
                                                "GO:0006536: glutamate metabolic process" ,
                                                "GO:0046890: regulation of lipid biosynthetic process",
                                                "GO:0006882: intracellular zinc ion homeostasis",
                                                "GO:0010043: response to zinc ion",
                                                "GO:0046688: response to copper ion",
                                                "GO:0009725: response to hormone",
                                                "GO:0071276: cellular response to cadmium ion",
                                                "GO:0007167: enzyme-linked receptor protein signaling pathway",
                                                "GO:0032922: circadian regulation of gene expression",
                                                "GO:0019430: removal of superoxide radicals")))

ggplot(GO) +
  geom_point(aes(x = group, 
                 y = Process,
                 size = (1-pval),
                 color = Fold_enrichment)) +
  scale_size(limits=c(.95,1),
             breaks=c(.95, .96,.97,.98,.99,1),
             labels=c(".05",".04",".03",".02",".01","<.01"),
             name = "p-value",
             guide="legend")+
  scale_colour_viridis(option = 'plasma', end = 0.85, begin = 0.2, limits = c(0,50), oob = squish,
                       labels = c('0', '10', '20', '30', '40', '50+'), name = 'Fold enrichment') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  #scale_x_discrete(drop = FALSE)+
  scale_y_discrete(limits=rev)

## Dot plot just stromal lineage ----
GO <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_TreatmentWithinCluster_GOresults.tsv')
GO$Treatment <- sub('.*\\_', '', GO$group)
GO$Cluster <- sub("_[^_]+$", "", GO$group)
vec <- c('Endothelial_35', # Endothelial
         'Fibroblast_31') 
GO <- GO[GO$Cluster %in% vec,] 
GO$pval <- as.numeric(GO$pval)
GO$Fold_enrichment <- as.numeric(GO$Fold_enrichment)
GO <- GO[GO$GO.ID %in% c(sGO_PEDV$GO.ID),]
table(GO$GO.ID, GO$Treatment)
GO$Process <- paste(GO$GO.ID, GO$Term, sep = ': ')

GO$group <- factor(GO$group, levels = c(paste(c('Endothelial_35', # Endothelial
                                                'Fibroblast_31'
), 'PEDV', sep = '_'),
paste(c('Endothelial_35', # Endothelial
        'Fibroblast_31'
), 'Mock', sep = '_')))
GO$Process <- factor(GO$Process, levels = rev(c("GO:0046718: viral entry into host cell",
                                            "GO:0002474: antigen processing and presentation of peptide antigen via MHC class I",
                                            "GO:0062208: positive regulation of pattern recognition receptor signaling pathway",
                                            "GO:0002753: cytosolic pattern recognition receptor signaling pathway",
                                            "GO:0071357: cellular response to type I interferon",
                                            "GO:0032727: positive regulation of interferon-alpha production",
                                            "GO:0032728: positive regulation of interferon-beta production",
                                            "GO:0060333: type II interferon-mediated signaling pathway",
                                            "GO:0031398: positive regulation of protein ubiquitination",
                                            "GO:0070936: protein K48-linked ubiquitination",
                                            "GO:2000147: positive regulation of cell motility",
                                            "GO:0030335: positive regulation of cell migration")))

ggplot(GO) +
  geom_point(aes(x = group, 
                 y = Process,
                 size = (1-pval),
                 color = Fold_enrichment)) +
  scale_size(limits=c(.95,1),
             breaks=c(.95, .96,.97,.98,.99,1),
             labels=c(".05",".04",".03",".02",".01","<.01"),
             name = "p-value",
             guide="legend")+
  scale_colour_viridis(option = 'plasma', end = 0.85, begin = 0.2, limits = c(0,50), oob = squish,
                       labels = c('0', '10', '20', '30', '40', '50+'), name = 'Fold enrichment') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  #scale_x_discrete(drop = FALSE)+
  scale_y_discrete(limits=rev)

# Session info ----
report(sessionInfo())
