library(Seurat)
library(SeuratDisk)
library(dplyr)
library(readxl)
library(writexl)
library(report)

# Load Seurat object ----
seu <- LoadH5Seurat('/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples.h5Seurat')

# Add meta data ----

seu$ID_JPP <- seu$orig.ident
Idents(seu) <- seu$ID_JPP
seu <- RenameIdents(object = seu, 
                    '1' = 'JPP1',
                    '2' = 'JPP3',
                    '3' = 'JPP4',
                    '4' = 'JPP5',
                    '5' = 'JPP6',
                    '6' = 'JPP7',
                    '7' = 'JPP8',
                    '8' = 'JPP9')
seu$ID_JPP <- Idents(seu) 

seu$AnimalID <- seu$orig.ident
Idents(seu) <- seu$AnimalID
seu <- RenameIdents(object = seu, 
                    '1' = '926',
                    '2' = '928',
                    '3' = '929',
                    '4' = '921',
                    '5' = '922',
                    '6' = '923',
                    '7' = '924',
                    '8' = '925')
seu$AnimalID <- Idents(seu)

seu$Treatment <- seu$orig.ident
Idents(seu) <- seu$Treatment
seu <- RenameIdents(object = seu, 
                    '1' = 'Mock',
                    '2' = 'Mock',
                    '3' = 'Mock',
                    '4' = 'PEDV',
                    '5' = 'PEDV',
                    '6' = 'PEDV',
                    '7' = 'PEDV',
                    '8' = 'PEDV')
seu$Treatment <- Idents(seu)

seu$LocalInfection_IHC <- seu$orig.ident
Idents(seu) <- seu$LocalInfection_IHC
seu <- RenameIdents(object = seu, 
                    '1' = 'Uninfected',
                    '2' = 'Uninfected',
                    '3' = 'Uninfected',
                    '4' = 'PEDV',
                    '5' = 'PEDV',
                    '6' = 'PEDV',
                    '7' = 'Uninfected',
                    '8' = 'PEDV')
seu$LocalInfection_IHC <- Idents(seu)

seu$JPPDistance_cm <- seu$orig.ident
Idents(seu) <- seu$JPPDistance_cm
seu <- RenameIdents(object = seu, 
                    '1' = 480,
                    '2' = 100,
                    '3' = 280,
                    '4' = 100,
                    '5' = 380,
                    '6' = 320,
                    '7' = 150,
                    '8' = 270)
seu$JPPDistance_cm <- Idents(seu)

seu$ClinicalScore_0dpi <- seu$orig.ident
Idents(seu) <- seu$ClinicalScore_0dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 0,
                    '5' = 0,
                    '6' = 0,
                    '7' = 0,
                    '8' = 0)
seu$ClinicalScore_0dpi <- Idents(seu)

seu$ClinicalScore_1dpi <- seu$orig.ident
Idents(seu) <- seu$ClinicalScore_1dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 1,
                    '5' = 0,
                    '6' = 0,
                    '7' = 0,
                    '8' = 1)
seu$ClinicalScore_1dpi <- Idents(seu)

seu$ClinicalScore_2dpi <- seu$orig.ident
Idents(seu) <- seu$ClinicalScore_2dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 1,
                    '5' = 2,
                    '6' = 0,
                    '7' = 0,
                    '8' = 0)
seu$ClinicalScore_2dpi <- Idents(seu)

seu$ClinicalScore_3dpi <- seu$orig.ident
Idents(seu) <- seu$ClinicalScore_3dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 1,
                    '5' = 2,
                    '6' = 1,
                    '7' = 0,
                    '8' = 2)
seu$ClinicalScore_3dpi <- Idents(seu)

seu$ClinicalScore_4dpi <- seu$orig.ident
Idents(seu) <- seu$ClinicalScore_4dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 1,
                    '5' = 2,
                    '6' = 1,
                    '7' = 0,
                    '8' = 1)
seu$ClinicalScore_4dpi <- Idents(seu)

seu$ClinicalScore_5dpi <- seu$orig.ident
Idents(seu) <- seu$ClinicalScore_5dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 1,
                    '5' = 2,
                    '6' = 1,
                    '7' = 0,
                    '8' = 2)
seu$ClinicalScore_5dpi <- Idents(seu)

seu$RectalSwabViralCopies_0dpi <- seu$orig.ident
Idents(seu) <- seu$RectalSwabViralCopies_0dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 0,
                    '5' = 0,
                    '6' = 0,
                    '7' = 0,
                    '8' = 0)
seu$RectalSwabViralCopies_0dpi <- Idents(seu)

seu$RectalSwabViralCopies_1dpi <- seu$orig.ident
Idents(seu) <- seu$RectalSwabViralCopies_1dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 0,
                    '5' = 21.76076,
                    '6' = 0,
                    '7' = 0,
                    '8' = 0)
seu$RectalSwabViralCopies_1dpi <- Idents(seu)

seu$RectalSwabViralCopies_2dpi <- seu$orig.ident
Idents(seu) <- seu$RectalSwabViralCopies_2dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 66.09873199,
                    '5' = 76156.01563,
                    '6' = 54486.17969,
                    '7' = 26.23313904,
                    '8' = 8.501250267)
seu$RectalSwabViralCopies_2dpi <- Idents(seu)

seu$RectalSwabViralCopies_3dpi <- seu$orig.ident
Idents(seu) <- seu$RectalSwabViralCopies_3dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 950.2220459,
                    '5' = 91111.20313,
                    '6' = 436119.8125,
                    '7' = 609484.25,
                    '8' = 933.2930298)
seu$RectalSwabViralCopies_3dpi <- Idents(seu)

seu$RectalSwabViralCopies_4dpi <- seu$orig.ident
Idents(seu) <- seu$RectalSwabViralCopies_4dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 881605.5,
                    '5' = 169836.2344,
                    '6' = 30580.23047,
                    '7' = 285.620636,
                    '8' = 767.1186523)
seu$RectalSwabViralCopies_4dpi <- Idents(seu)

seu$RectalSwabViralCopies_5dpi <- seu$orig.ident
Idents(seu) <- seu$RectalSwabViralCopies_5dpi
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 931922.125,
                    '5' = 70322.11719,
                    '6' = 152286.625,
                    '7' = 441554.875,
                    '8' = 1436.541016)
seu$RectalSwabViralCopies_5dpi <- Idents(seu)

seu$VillusBlunting <- seu$orig.ident
Idents(seu) <- seu$VillusBlunting
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 3,
                    '5' = 4,
                    '6' = 3,
                    '7' = 0,
                    '8' = 1)
seu$VillusBlunting <- Idents(seu)

seu$VillusFusion <- seu$orig.ident
Idents(seu) <- seu$VillusFusion
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 3,
                    '5' = 4,
                    '6' = 3,
                    '7' = 0,
                    '8' = 2)
seu$VillusFusion <- Idents(seu)

seu$CryptElongation <- seu$orig.ident
Idents(seu) <- seu$CryptElongation
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 3,
                    '5' = 4,
                    '6' = 3,
                    '7' = 0,
                    '8' = 2)
seu$CryptElongation <- Idents(seu)

seu$EpithelialAttenuation <- seu$orig.ident
Idents(seu) <- seu$EpithelialAttenuation
seu <- RenameIdents(object = seu, 
                    '1' = 0,
                    '2' = 0,
                    '3' = 0,
                    '4' = 0,
                    '5' = 0,
                    '6' = 0,
                    '7' = 0,
                    '8' = 0)
seu$EpithelialAttenuation <- Idents(seu)

# Save object ----

SaveH5Seurat(seu, 
             filename = "/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/NormIntRedOutputs/AllSamples.h5Seurat", overwrite = TRUE)


### View session information ----
report(sessionInfo())
