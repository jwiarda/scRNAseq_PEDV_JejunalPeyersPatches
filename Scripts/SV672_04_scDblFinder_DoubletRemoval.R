library(scDblFinder)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(report)

## Perform doublet calculation within each sample
# Recommendation: best to perform this on data that has had empty drops removed but has not been filtered any further: https://github.com/plger/scDblFinder

### First make sure we have enough read counts in our cells to accurately detect doublets:
# Create Seurat object ----
data_dir <- c('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP1strainedCounts',
              '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP3strainedCounts',
              '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP4strainedCounts',
              '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP5strainedCounts',
              '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP6strainedCounts',
              '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP7strainedCounts',
              '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP8strainedCounts',
              '/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP9strainedCounts')
scRNA_data <- Read10X(data.dir = data_dir) # read the 10X data from all samples into a data matrix
seu = CreateSeuratObject(counts = scRNA_data)

# Plot genes/cell and UMIs/cell detected in each sample ----
meta <- seu@meta.data
ggplot(meta, aes(x=nCount_RNA,y=..density..)) + 
  geom_histogram(fill="white",color="black",bins=500) + 
  scale_x_continuous(breaks = seq(0, 50000, 2000), lim = c(0, 50000)) + # lots of data goes past our set upper limit
  facet_wrap(~orig.ident) +
  geom_vline(aes(xintercept=500),color="red",lty="longdash") + # move this cutoff line where you see fit
  RotatedAxis() + 
  ggtitle('nCount_RNA')
## From the above plot, we determine most of our cells have sufficiently high total transcript counts to perform doublet calculations, but we do have a handful of low-read cells.
## We opt to get rid of cells with <500 total UMIs before performing doublet removal
# Note: we should not be doing stringent cell QC filtering yet, JUST REMOVING EXTREMELY LOW READ DROPLETS. 
  # Since we only red in droplets with at least 500 UMIs from CellRanger to SoupX, we shouldn't have very many extremely low UMI cells now unless SoupX did a ton of ambient RNA removal on some droplets

# Filter out exceptionally low-transcript cells ----
select <- WhichCells(seu, expression = nCount_RNA > 500)
seu <- subset(seu, cells = select)

### Run doublet calculation on each sample individually ----
# Convert filtered Seurat object to SingleCellExperiment object:
sce <- as.SingleCellExperiment(seu)

# Calculate doublets:
set.seed(123) # set seed for reproducibility
sce <- scDblFinder(sce, 
                   samples="orig.ident", #run each sample individually
                   #BPPARAM=MulticoreParam(3), 
                   clusters=TRUE) # use clustering approach; clusters will be generated here since they have not been pre-calculated; clustering approach suggested for well-segmented data... non-clustering recommended for data with poor segregation, like a developmental trajectory: https://github.com/plger/scDblFinder
table(sce$scDblFinder.class, sce$orig.ident) # see table of predicted doublets
#1    2    3    4    5    6    7    8
#singlet 1593 8236 6986 7435 7516 9339 8238 4784
#doublet   88  860  802  795 1065 1361  837  472

# Convert SingleCellExperiment back to a Seurat object that includes doublet information:
logcounts(sce) <- assay(sce, "counts") # not actually log counts. See workaround at https://github.com/satijalab/seurat/issues/3746
seu <- as.Seurat(sce)
table(seu$scDblFinder.class)

### Save data ----
# Save Seurat object containing doublet information:
dir.create('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/scDblFinderOutputs/')
SaveH5Seurat(seu, 
             filename = "/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/scDblFinderOutputs/AllSamples_DoubletsCalculated.h5Seurat", overwrite = TRUE)

### View session information ----
report(sessionInfo())