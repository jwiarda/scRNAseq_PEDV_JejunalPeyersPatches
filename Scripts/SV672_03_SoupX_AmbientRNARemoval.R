library(DropletUtils) 
library(ggplot2) 
library(SoupX) 
library(report)

## Perform ambient RNA estimation and removal on each individual sample 

## Create a directory for results
dir.create('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/')

### Sample JPP1 ----
#### Provide file path to Cell Ranger outputs:
sc = load10X('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts/JPP1/outs') 

#### Check the meta data to make sure we also have cluster information and t-SNE coordinates:
head(sc$metaData, n = 3)

#### Plot data according to Cell Ranger cluster assignments and t-SNE coordinates:
dd = sc$metaData # create an object with all the metadata
mids = aggregate(cbind(tSNE1,tSNE2) ~ clusters,data=dd,FUN=mean) # determine t-SNE coordinates for middle of each cluster
gg = ggplot(dd,aes(tSNE1,tSNE2)) + # make a t-SNE plot
  geom_point(aes(colour=factor(clusters)),size=0.2) +
  geom_label(data=mids,aes(label=clusters)) 
plot(gg) # show plot

#### Check expression patterns for some canonical genes:
dd$CD3E = sc$toc["CD3E", ] # make column of gene expression values for CD3E (T cell gene)
dd$IgLambdaV = sc$toc["ENSSSCG00000038719", ] # make column of gene expression values for gene that codes for Ig lambda V region (B cell gene)
dd$CD79B = sc$toc["CD79B", ] # make column of gene expression values for CD79B (B cell gene)
dd$FABP6 = sc$toc["FABP6", ] # make column of gene expression values for gene that codes for FABP6 (epithelial cell gene)
dd$EPCAM = sc$toc["EPCAM", ] # make column of gene expression values for gene that codes for EPCAM (epithelial cell gene)
dd$GNLY = sc$toc["GNLY", ] # make column of gene expression values for gene that codes for GNLY (cytotoxicty gene)
dd$HBB = sc$toc["HBB", ] # make column of gene expression values for gene that codes for HBB (erythrocyte gene)

a1 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD3E > 0)) + ggtitle('CD3E') # which cells express this gene?
a2 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = IgLambdaV > 0)) + ggtitle('IgLambdaV')
a3 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = CD79B > 0)) + ggtitle('CD79B')
a4 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = FABP6 > 0)) + ggtitle('FABP6')
a5 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = EPCAM > 0)) + ggtitle('EPCAM')
a6 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = GNLY > 0)) + ggtitle('GNLY')
a7 <- ggplot(dd, aes(tSNE1,tSNE2)) + geom_point(aes(colour = HBB > 0)) + ggtitle('HBB')
(a1+a2+a3) / (a4+a5+a6) # show 6 plots

a1 <- plotMarkerMap(sc, "CD3E") + ggtitle('CD3E')  # if we assumed all cells were nothing but soup, which cells still show higher than expected expression for the gene (TRUE = expression levels higher than expected if cell was just soup, so likely real expression). This just gives us an idea of soup expression, this is NOT a formal analysis used for removing the soup RNA.
a2 <- plotMarkerMap(sc, "ENSSSCG00000038719") + ggtitle('IgLambdaV')
a3 <- plotMarkerMap(sc, "CD79B") + ggtitle('CD79B')
a4 <- plotMarkerMap(sc, "FABP6") + ggtitle('FABP6')
a5 <- plotMarkerMap(sc, "EPCAM") + ggtitle('EPCAM')
a6 <- plotMarkerMap(sc, "GNLY") + ggtitle('GNLY')
a7 <- plotMarkerMap(sc, "HBB") + ggtitle('HBB')
(a1+a2+a3) / (a4+a5+a6) # show 6 plots

#What we see in these plots is some misplaced gene expression, indicating we have RNA soup to remove!

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### See how soup removal affects the genes we assessed expression patterns for earlier:
a1 <- plotChangeMap(sc, out, "CD3E") + ggtitle('CD3E')
a2 <- plotChangeMap(sc, out, "ENSSSCG00000038719") + ggtitle('IgLambdaV')
a3 <- plotChangeMap(sc, out, "CD79B") + ggtitle('CD79B')
a4 <- plotChangeMap(sc, out, "FABP6") + ggtitle('FABP6')
a5 <- plotChangeMap(sc, out, "EPCAM") + ggtitle('EPCAM')
a6 <- plotChangeMap(sc, out, "GNLY") + ggtitle('GNLY')
a7 <- plotChangeMap(sc, out, "HBB") + ggtitle('HBB')
(a1+a2+a3) / (a4+a5+a6) # show 6 plots

#### Save our strained count matrix to a new location:
write10xCounts("/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP1strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

#Now I do the same for the rest of the files, but I won't show most of the output plots to save some space (though I did look at these plots as a double check)

### Sample JPP3 ----
#### Provide file path to Cell Ranger outputs:
sc = load10X('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts/JPP3/outs') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP3strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample JPP4 ----
#### Provide file path to Cell Ranger outputs:
sc = load10X('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts/JPP4/outs') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP4strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample JPP5 ----
#### Provide file path to Cell Ranger outputs:
sc = load10X('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts/JPP5/outs') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP5strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample JPP6 ----
#### Provide file path to Cell Ranger outputs:
sc = load10X('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts/JPP6/outs') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP6strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample JPP7 ----
#### Provide file path to Cell Ranger outputs:
sc = load10X('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts/JPP7/outs') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP7strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample JPP8 ----
#### Provide file path to Cell Ranger outputs:
sc = load10X('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts/JPP8/outs') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP8strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### Sample JPP9 ----
#### Provide file path to Cell Ranger outputs:
sc = load10X('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts/JPP9/outs') 

#### Calculate the RNA soup fraction:
sc = autoEstCont(sc) # estimate the fraction of RNAs belonging to soup
out = adjustCounts(sc) # create a corrected count matrix

#### See which genes were most affected by our soup correction:
cntSoggy = rowSums(sc$toc > 0) # list cells with counts greater than 0 before correction for each gene
cntStrained = rowSums(out > 0) # list cells with counts greater than 0 after correction for each gene
tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10) # list the 10 most affected genes that had expression reduced in total # of cells
tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 10) # list the 10 genes that had greatest overall quantities reduced

#### Save our strained count matrix to a new location:
write10xCounts("/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/SoupX/JPP9strainedCounts", out, version = "3", overwrite = TRUE)
rm(dd, gg, mids, out, sc, cntSoggy, cntStrained)

### View session information ----
report(sessionInfo())
