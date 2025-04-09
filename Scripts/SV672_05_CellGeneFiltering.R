library(DropletUtils) 
library(Seurat)
library(readxl)
library(ggplot2)
library(scales)
library(SeuratDisk)
library(report)

# Load Seurat object of all samples with doublets calculated ----
## note that we have already filtered out some cells that had <500 UMIs per cell....these were filtered out BEFORE doublet calculation
## At this point, doublets have NOT been removed, only calculated
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/scDblFinderOutputs/AllSamples_DoubletsCalculated.h5Seurat')

# See how many cells/genes in Seurat object:
seu

#See how many cells in each sample:
table(seu$orig.ident)

#We need to calculate the percentage of mitochondrial genes expressed within each cell. ----

##To start, we need to form a list of mitochondrial genes to look at:
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

### Plot QC metrics ----
#Violin plots:
VlnPlot(seu, 
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), # QC metric metadata to plot
        group.by = 'orig.ident', # plot each sample separately
        pt.size = 0.01, 
        ncol = 3)

#Histograms with our desired thresholds shown from looking at all plots:
meta <- seu@meta.data
ggplot(meta, aes(x=percent_mito,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 50, 5), lim = c(0, 50)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=15),color="red",lty="longdash") + # move this cutoff line where you see fit
  RotatedAxis() + 
  ggtitle('percent_mito')

ggplot(meta, aes(x=nFeature_RNA,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 4000, 250), lim = c(0, 4000)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=250),color="red",lty="longdash") + # move this cutoff line where you see fit
  RotatedAxis() + 
  ggtitle('nFeature_RNA')

ggplot(meta, aes(x=nCount_RNA,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 60000, 2000), lim = c(0, 60000)) + 
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=3000),color="red",lty="longdash") + # move this cutoff line where you see fit
  RotatedAxis() + 
  ggtitle('nCount_RNA')

# Remake some xy plots with cutoffs:
#Number of genes vs. percent mitochondrial reads (color = # UMIs):
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA, color = nCount_RNA))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() +
  scale_colour_gradient(low = "gold", high = "red", limits=c(0, 4000), oob=squish) + 
  geom_hline(aes(yintercept=250),color="blue",lty="longdash") + 
  geom_vline(aes(xintercept=15),color="blue",lty="longdash")

#Number of genes vs. percent mitochondrial reads (color = doublet classification):
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA, color = scDblFinder.class))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  geom_hline(aes(yintercept=250),color="blue",lty="longdash") + 
  geom_vline(aes(xintercept=15),color="blue",lty="longdash")

# Number of genes vs. number reads (color = mitochondrial reads):
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA, color = percent_mito))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,60000) +
  ylim(0,8000) +
  scale_colour_gradient(low = "gold", high = "red", limits=c(0, 20), oob=squish) + 
  geom_hline(aes(yintercept=250),color="blue",lty="longdash") + 
  geom_vline(aes(xintercept=3000),color="blue",lty="longdash")

# Number of genes vs. number reads (color = doublet classification):
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA, color = scDblFinder.class))+
  geom_point(size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,60000) +
  ylim(0,8000) + 
  geom_hline(aes(yintercept=250),color="blue",lty="longdash") + 
  geom_vline(aes(xintercept=3000),color="blue",lty="longdash")


#Established QC thresholds are:
##Keep cells with <15% mitochondrial reads
##Keep cells with >250 genes detected
##Keep cells with >1500 total reads
##Keep cells not called as doublets

##Filter out poor quality cells ----
#Identify cells passing each/every QC filter:
keepMito <- WhichCells(seu, expression = percent_mito < 15) # cells passing mito filter
length(keepMito) # how many cells pass this QC filter?

keepGenes <- WhichCells(seu, expression = nFeature_RNA > 250) # cells passing gene filter
length(keepGenes) # how many cells pass this QC filter?

keepUMI <- WhichCells(seu, expression = nCount_RNA > 3000) # cells passing UMI filter
length(keepUMI) # how many cells pass this QC filter?

keepDub <- WhichCells(seu, expression = scDblFinder.class == 'singlet') # cells passing doublet filter
length(keepDub) # how many cells pass this QC filter?

keep <- Reduce(intersect, list(keepMito, keepGenes, keepUMI, keepDub)) # cells passing all QC filters
length(keep) # how many cells pass this QC filter?

rm(keepMito, keepGenes, keepUMI, keepDub) # free up space

# Create new Seurat object with only cells passing all QC filters: ----
seuKeep <- subset(seu, 
                  cells = keep)

# How many cells in each sample now?
table(seuKeep$orig.ident)
#  1    2    3    4    5    6    7    8 
# 1295 6647 5462 6310 5967 7656 6967 3947 

### Look over QC plots to see which cells did or did not pass QC ----
# Number of genes vs. percent mitochondrial reads:
metaKeep <- seuKeep@meta.data
ggplot(meta, aes(x=percent_mito, y=nFeature_RNA))+
  geom_point(color = 'black', size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() +
  geom_point(data=metaKeep, aes(x=percent_mito, y=nFeature_RNA), color = 'red', size = 0.75)

# Number of genes vs. number reads:
ggplot(meta, aes(x=nCount_RNA, y=nFeature_RNA))+
  geom_point(color = 'black', size = 0.75) + 
  facet_wrap(~orig.ident, nrow =1)+
  theme_get() + 
  xlim(0,60000) +
  ylim(0,8000) +
  geom_point(data=metaKeep, aes(x=nCount_RNA, y=nFeature_RNA), color = 'red', size = 0.75)

## Save counts of cells passing QC from each sample ----
counts <- seuKeep@assays$RNA@counts
seuNew <- CreateSeuratObject(counts = counts, 
                             min.cells = 1) # include only genes expressed in at least one cell; THIS SUFFICES AS THE GENE FILTERING STEP FOR OUR DATA
seuNew # see how many genes/cells in new Seurat object

# Save Seurat object containing filtered cells: ----
seuNew[["RNA"]] <- as(object = seuNew[["RNA"]], Class = "Assay") # workaround at https://github.com/mojaveazure/seurat-disk/issues/27
dir.create('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/FilterQCOutputs/')
SaveH5Seurat(seuNew, 
             filename = "/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/FilterQCOutputs/AllSamples_Filtered.h5Seurat", overwrite = TRUE)

### View session information ----
report(sessionInfo())
