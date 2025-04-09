#!/bin/bash

#SBATCH --job-name=CellRanger                          # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 3                                            # number of nodes in this job
#SBATCH -n 40                                             # number of cores/tasks in this job
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=200G   					# memory e.g.: 100G ; 250G ; 100M etc... (need mem node for some reqs)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jayne.wiarda@usda.gov					# enter your email address here USDA emails sometimes dont work, I use gmail
#SBATCH -A nadc_prrsv

# ENTER COMMANDS HERE:

# First install CellRanger according to instructions on 10X Genomics website.
# CellRanger v8.0.1

# Define CellRanger path to use for correct version:
export PATH=/home/jayne.wiarda/yard/apps/cellranger-8.0.1:$PATH

# Perform alignment and counting to the reference genome for each sample
mkdir /project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts
cd /project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerCounts

cellranger count --id=JPP1 \
                   --transcriptome=/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex/Sscrofa11-1 \
                   --fastqs=/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane1,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane2 \
                   --sample=JPP1 \
                   --create-bam=true \
                   --localcores=38

#JPP3 sample had some deeper sequencing done on another flow cell run, so we have more lanes to read in for this particular sample
cellranger count --id=JPP3 \
                   --transcriptome=/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex/Sscrofa11-1 \
                   --fastqs=/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane1,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane2,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_241125/lane1,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_241125/lane2 \
                   --sample=JPP3 \
                   --create-bam=true \
                   --localcores=38

cellranger count --id=JPP4 \
                   --transcriptome=/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex/Sscrofa11-1 \
                   --fastqs=/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane1,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane2 \
                   --sample=JPP4 \
                   --create-bam=true \
                   --localcores=38

cellranger count --id=JPP5 \
                   --transcriptome=/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex/Sscrofa11-1 \
                   --fastqs=/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane1,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane2 \
                   --sample=JPP5 \
                   --create-bam=true \
                   --localcores=38

cellranger count --id=JPP6 \
                   --transcriptome=/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex/Sscrofa11-1 \
                   --fastqs=/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane1,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane2 \
                   --sample=JPP6 \
                   --create-bam=true \
                   --localcores=38
                   
cellranger count --id=JPP7 \
                   --transcriptome=/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex/Sscrofa11-1 \
                   --fastqs=/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane1,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane2 \
                   --sample=JPP7 \
                   --create-bam=true \
                   --localcores=38
                   
cellranger count --id=JPP8 \
                   --transcriptome=/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex/Sscrofa11-1 \
                   --fastqs=/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane1,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane2 \
                   --sample=JPP8 \
                   --create-bam=true \
                   --localcores=38
                   
cellranger count --id=JPP9 \
                   --transcriptome=/project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex/Sscrofa11-1 \
                   --fastqs=/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane1,/project/nadc_prrsv/Wiarda/SequencingFiles/wiarda_240717/lane2 \
                   --sample=JPP9 \
                   --create-bam=true \
                   --localcores=38
                   
#End of file
