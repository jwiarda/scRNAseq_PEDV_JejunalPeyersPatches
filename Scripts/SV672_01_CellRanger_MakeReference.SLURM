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

# Add PEDV genome to annotation file:

# FASTA from https://ncbi.nlm.nih.gov/nuccore/KM975737.1?report=fasta; make sure to save as a simple text file, not a rich text file (RTF)
cd /project/nadc_prrsv/Wiarda/GeneAnnotationFiles
# ran the following command to find out there are 28038 bases in the file, which will be needed as an input value later
#cat PEDV_USA_NC_2013_49469.fa  | grep -v "^>" | tr -d "\n" | wc -c
# make custom GTF file for gene:
echo -e 'PEDV\tunknown\texon\t1\t28038\t.\t+\t.\tgene_id "PEDV"; transcript_id "PEDV"; gene_name "PEDV"; gene_biotype "protein_coding";' > PEDV_USA_NC_2013_49469.gtf
cp Sus_scrofa.Sscrofa11.1.dna.toplevel.fa Sus_scrofa.Sscrofa11.1.dna.toplevel_PEDV_USA_NC_2013_49469.fa
cat PEDV_USA_NC_2013_49469.fa >> Sus_scrofa.Sscrofa11.1.dna.toplevel_PEDV_USA_NC_2013_49469.fa
cp Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.gtf Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS_PEDV_USA_NC_2013_49469.gtf 
cat PEDV_USA_NC_2013_49469.gtf  >> Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS_PEDV_USA_NC_2013_49469.gtf 


# First install CellRanger according to instructions on 10X Genomics website.
# CellRanger v8.0.1

# Define CellRanger path to use for correct version:
export PATH=/home/jayne.wiarda/yard/apps/cellranger-8.0.1:$PATH

# Make reference genome:
# Obtained genome file http://ftp.ensembl.org/pub/release-97/fasta/sus_scrofa/dna/ & annotation file modified from http://ftp.ensembl.org/pub/release-97/gtf/sus_scrofa/

mkdir /project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex
cd /project/nadc_prrsv/Wiarda/SV672_scRNAseq_PEDVJejPP/CellRangerGenomeIndex

cellranger mkref --genome=Sscrofa11-1 \
                 --fasta=/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/Sus_scrofa.Sscrofa11.1.dna.toplevel_PEDV_USA_NC_2013_49469.fa \
                 --genes=/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS_PEDV_USA_NC_2013_49469.gtf \
                 --ref-version=11.1.97 \
		 --nthreads=8

#End of file
