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
# Load Seurat object ----
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/trajectorySubset_PolarizedEnterocytes.h5seurat')
DefaultAssay(seu) <- 'SCT'

## Read in & filter DE data ----
AllCellsDGE <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/DGE/DEGs_ClustersPairwise_PEDVEnterocyteInfectionStates_unflitered.tsv')
AllCellsDGE <- subset(AllCellsDGE, logFC > 0.25 & padj < 0.05)
AllCellsDGE <- subset(AllCellsDGE, pct_in > 0.1 | pct_out > 0.1)

Hs_BystanderVsMock <- subset(pwAll, group == 'HomeostaticEnterocyte_Bystander' & (pop1 == 'HomeostaticEnterocyte_Mock' | pop2 == 'HomeostaticEnterocyte_Mock'))
Hs_BystanderVsMock$group <- 'Hs_BystanderVsMock'
Hs_InfectedVsMock <- subset(pwAll, group == 'HomeostaticEnterocyte_Infected' & (pop1 == 'HomeostaticEnterocyte_Mock' | pop2 == 'HomeostaticEnterocyte_Mock'))
Hs_InfectedVsMock$group <- 'Hs_InfectedVsMock'

St_BystanderVsMock <- subset(pwAll, group == 'StressedEnterocyte_Bystander' & (pop1 == 'StressedEnterocyte_Mock' | pop2 == 'StressedEnterocyte_Mock'))
St_BystanderVsMock$group <- 'St_BystanderVsMock'
St_InfectedVsMock <- subset(pwAll, group == 'StressedEnterocyte_Infected' & (pop1 == 'StressedEnterocyte_Mock' | pop2 == 'StressedEnterocyte_Mock'))
St_InfectedVsMock$group <- 'St_InfectedVsMock'

# Conserved signature ----
ConservedPEDVSignature <- intersect(St_InfectedVsMock$feature, (intersect(St_BystanderVsMock$feature, (intersect(Hs_BystanderVsMock$feature, Hs_InfectedVsMock$feature)))))
length(ConservedPEDVSignature)
AllCellsDGE <- data.frame(ConservedPEDVSignature)
colnames(AllCellsDGE) <- 'feature'
AllCellsDGE$group <- 'ConservedPEDVSignature'

# GO analysis ----
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
    
    all_GO <- read_tsv('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_gene_to_GO_enterocytes.tsv')
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
  write_tsv('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_gene_to_GO_enterocytes.tsv')

# Cell type-specific GO terms ----

AllCells_results <- 
  AllCellsDGE %>% 
  #filter(logFC > 0) %>% 
  group_by(group) %>% 
  nest() %>% 
  mutate(enriched_genes=map(data, ~.x %>% pull(feature))) %>% 
  mutate(GO_results=
           map(enriched_genes, ~topGO_wrapper(myInterestingGenes = .x,
                                              mapping_file = '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_gene_to_GO_enterocytes.tsv'))) %>% 
  mutate(filt_results=purrr::map(.x = GO_results, .f = ~filter(.x, pval < 0.05))) %>% 
  dplyr::select(group, filt_results) %>% 
  unnest(cols=filt_results)

AllCells_results$Fold_enrichment <- AllCells_results$Significant/AllCells_results$Expected # add in stat for fold enrichment
AllCells_results <- subset(AllCells_results, Significant > 1) # remove GO terms with only one enriched genedim

AllCells_results %>% write_tsv('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_EnterocyteInfectionState_ConservedPEDVSignature_GOresults.tsv')

# Identify top processes ----
# Criteria used was GO term had pval < 0.05 in at least half the populations of a lineage, then took the top terms with lowest p value. 

## Plot GO terms ----
GO <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_EnterocyteInfectionState_ConservedPEDVSignature_GOresults.tsv')
GO$pval <- as.numeric(GO$pval)

ag <- slice_min(GO, GO$pval, n = 20)
newGO <- GO[GO$GO.ID %in% ag$GO.ID,]

# Create dot plot ----
GO <- read.delim('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/GeneOntology/AllCells_EnterocyteInfectionState_ConservedPEDVSignature_GOresults.tsv')
GO <- GO[GO$GO.ID %in% c(newGO$GO.ID),]
table(GO$GO.ID, GO$group)
GO$Process <- paste(GO$GO.ID, GO$Term, sep = ': ')
GO$Process <- factor(GO$Process, levels = rev(c('GO:0071357: cellular response to type I interferon',
                                                'GO:0045087: innate immune response',
                                                'GO:0062208: positive regulation of pattern recognition receptor signaling pathway',
                                                'GO:0035455: response to interferon-alpha',
                                                'GO:0039535: regulation of RIG-I signaling pathway',
                                                'GO:0002888: positive regulation of myeloid leukocyte mediated immunity',
                                                'GO:0032728: positive regulation of interferon-beta production',
                                                'GO:0034314: Arp2/3 complex-mediated actin nucleation',
                                                'GO:0099638: endosome to plasma membrane protein transport',
                                                'GO:0042274: ribosomal small subunit biogenesis',
                                                'GO:0071356: cellular response to tumor necrosis factor',
                                                'GO:0006123: mitochondrial electron transport, cytochrome c to oxygen',
                                                'GO:0071622: regulation of granulocyte chemotaxis',
                                                'GO:0050829: defense response to Gram-negative bacterium',
                                                'GO:0071709: membrane assembly',
                                                'GO:0002690: positive regulation of leukocyte chemotaxis',
                                                'GO:0034656: nucleobase-containing small molecule catabolic process',
                                                'GO:0010758: regulation of macrophage chemotaxis',
                                                'GO:0001919: regulation of receptor recycling',
                                                'GO:0072529: pyrimidine-containing compound catabolic process')))
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
  scale_colour_viridis(option = 'plasma', end = 0.85, begin = 0.2, limits = c(0,20), oob = squish,
                       name = 'Fold enrichment') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1)) +
  #scale_x_discrete(drop = FALSE)+
  scale_y_discrete(limits=rev)

ggplot(GO) +
  geom_bar(stat = "identity", aes(x = Process, 
                                  y = Fold_enrichment,
                                  fill = pval)) +
  coord_flip() +
  theme_bw()+
  scale_fill_viridis(option = 'plasma', end = 0.85, begin = 0.2, limits = c(0, 0.002), oob = squish,
                     name = 'P value', direction = -1)

## Plot genes from top GO terms ----
#seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/Trajectory/trajectorySubset_PolarizedEnterocytes.h5seurat')
#DefaultAssay(seu) <- 'SCT'
seu$Infection <- NA
seu$Infection[seu$Treatment == 'PEDV' & seu[["SCT"]]@data['PEDV',] > 0] <- 'Infected'
seu$Infection[seu$Treatment == 'PEDV' & seu[["SCT"]]@data['PEDV',] == 0] <- 'Bystander'
seu$Infection[seu$Treatment == 'Mock'] <- 'Mock'

seu$pseudotime1[is.na(seu$pseudotime1)] <- -5
seu$pseudotime2[is.na(seu$pseudotime2)] <- -5
seu$EnterocytePolarization <- 'Both'
seu$EnterocytePolarization[seu$pseudotime1 > 0 & seu$pseudotime2 < 0] <- 'Traj1'
seu$EnterocytePolarization[seu$pseudotime2 > 0 & seu$pseudotime1 < 0] <- 'Traj2'

seu$EnterocyteState <- seu$EnterocytePolarization
seu$EnterocyteState[seu$EnterocytePolarization == 'Traj1'] <- 'HomeostaticEnterocyte'
seu$EnterocyteState[seu$EnterocytePolarization == 'Traj2'] <- 'StressedEnterocyte'
Idents(seu) <- seu$EnterocyteState
seu <- subset(seu, idents = c('HomeostaticEnterocyte', 'StressedEnterocyte'))

seu$EnterocyteInfection <- paste(seu$EnterocyteState, seu$Infection, sep = '_')
seu$EnterocyteInfection <- factor(seu$EnterocyteInfection, levels = rev(c('HomeostaticEnterocyte_Mock', 'StressedEnterocyte_Mock', 
                                                                          'HomeostaticEnterocyte_Bystander', 'StressedEnterocyte_Bystander', 
                                                                          'HomeostaticEnterocyte_Infected', 'StressedEnterocyte_Infected')))

#DotPlot(seu, group.by = 'EnterocyteInfection', features = ConservedPEDVSignature, cols = c('gold', 'darkmagenta')) + RotatedAxis()

DotPlot(seu, group.by = 'EnterocyteInfection', features = unique(c('ANXA2', 'ARPC1A', 'ARPC2', 'ARPC3', 'ARPC4', 'BST2', 'C3', 'CDA', 'CHMP5', 'COX4I1', 'COX5A', 'COX7A2', 'DDX58', 'DDX60', 'DHX58', 'DPP4', 'EIF2AK2', 'ENTPD7', 'ERBIN', 'F2RL1', 'FCN2', 
                                                                   'GAPDH', 'HERC5', 'IFI6', 'IFIT1', 'IFITM3', 'IL34', 'IRF7', 'ISG15', 'KRT8', 'LYPD8', 'MAPK3', 'MCU', 'MX1', 'NLRC5', 
                                                                   'OAS2', 'OASL', 'OPTN', 'PML', 'PTK2B', 'RAB11A', 'RPS12', 'RPS17', 'RPS3', 'RPS4X', 'RPS7', 'RPS8', 'RPSA', 'RSAD2', 'S100A10', 'SORL1', 'SPTBN1', 'STAT1', 'STAT2', 
                                                                   'TRIM22', 'TRIM40', 'UBE2L6', 'UPP1', 'USP15', 'USP18', 'USP19', 'ZDHHC18', 'ZNFX1')), cols = c('gold', 'darkmagenta')) + RotatedAxis() # just genes from top 10 GO processes

# Session info ----
report(sessionInfo())
