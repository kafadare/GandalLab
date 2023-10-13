##This is the code used to generate the adult genExp counts and genExp tpm, Sep 2023 EK
#Set up packages and functions
setwd("~/project-gandalm/GandalLab")
output_dir <- ("~/project-gandalm/comb_data/full_set/")
source("data_functions.R")
source("analysis_fcts.R")
# List of required packages
#required_packages <- c("magrittr", "genio", "dplyr", "BiocManager", "SummarizedExperiment")
required_packages <- c("magrittr","dplyr")
# Install or load missing packages
#load_install_pkg(required_packages)
library(magrittr)
library(dplyr)


adult_se <- readRDS("/u/home/k/kafadare/project-gandalm/AdultBigBrain_gene_exp_raw_042923.RDS") #filtered for lowly expressed genes
#adult_se <- readRDS("/u/project/gandalm/kafadare/AdultBigBrain_tx_exp_raw_042923.RDS") #transcripts, non-filtered
#gene_raw = tximeta::summarizeToGene(adult_se) #convert to gene counts from transcripts
#assays(adult_se)
adult <- list()
adult$genExp.counts <- adult_se %>% assay(.,1)
adult$genExp.tpm <-  adult_se %>% assay(.,2)
write.table(adult$raw_data.counts,file="adult.counts.scaled.tsv",quote=FALSE, sep='\t')
write.table(adult$raw_data.tpm,file="adult.TPM.tsv",quote=FALSE, sep='\t')