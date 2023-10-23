##Combine geneExp datasets and metadata for fetal and adult EUR

#Set up packages and functions
setwd("~/project-gandalm/GandalLab")
output_dir <- ("~/project-gandalm/comb_data/new_data/")
source("data_functions.R")
source("analysis_fcts.R")
# List of required packages
#required_packages <- c("magrittr", "genio", "dplyr", "BiocManager", "SummarizedExperiment")
required_packages <- c("magrittr","dplyr")
# Install or load missing packages
load_install_pkg(required_packages)
#library(magrittr)
#library(dplyr)
#library(SummarizedExperiment)
#library(argparser)
#library(tximport)
#library(tximeta)
library(SummarizedExperiment)


#fetal paths
fetal_counts_path <- "~/project-gandalm/fetal_rmlow_counts.tsv"
fetal_tpm_path <- "~/project-gandalm/fetal_rmlow_tpm.tsv"
#fetal_counts_path <- "/u/project/gandalm/cindywen/isoform_twas/salmon/expression.final/gene.noVersion.scaled.counts.tsv"
#fetal_tpm_path <- "/u/project/gandalm/cindywen/isoform_twas/salmon/expression.final/gene.noVersion.TPM.tsv"
fetal_meta_path <- "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/metadata_654.tsv"
#adult paths
adult_counts_path <- "~/project-gandalm/adult.counts.scaled.tsv"
adult_tpm_path <- "~/project-gandalm/adult.TPM.tsv"
adult_meta_path <- "/u/project/gandalm/kafadare/cov_hcp0_gene.txt"
#load data for fetal & adult
fetal <- load_data(names = c("counts", "tpm", "meta"), fetal_counts_path, fetal_tpm_path, fetal_meta_path)
adult <- load_data(names = c("counts", "tpm", "meta"), adult_counts_path, adult_tpm_path, adult_meta_path)
adult$ancestry <-  read.table("/u/project/gandalm/kafadare/pops.txt", header = F, sep = "\t", col.names = c("id", "ancestry"))
colnames(adult$counts) <- id_format_fix(colnames(adult$counts))
colnames(adult$tpm) <- id_format_fix(colnames(adult$tpm))

#fetal$tpm <- read.table("/u/project/gandalm/cindywen/isoform_twas/salmon/expression.final/gene.noVersion.TPM.tsv", header = TRUE, sep = "\t") #tpm file

#reformat adult metadata and add batch column from the RDS file
adult_se <- readRDS("/u/home/k/kafadare/project-gandalm/AdultBigBrain_gene_exp_raw_042923.RDS")
adult$meta <- as.data.frame(t(adult$meta))
colnames(adult$meta) <- adult$meta[1, ]
adult$meta <- adult$meta[-1, ]
#add adult data batch column in metadata df
batch_index <- match(colData(adult_se)$names,rownames(adult$meta))
adult$meta$study <- colData(adult_se)$batch[batch_index]

#get index match for ancestry data, append column to meta data
adult$ancestry$id <- id_format_fix(adult$ancestry$id)
missing_ids <- rownames(adult$meta)[!(rownames(adult$meta) %in% adult$ancestry$id)] # 61 ids "missing" from the ancestry file
adult$meta <- adult$meta[(rownames(adult$meta) %in% adult$ancestry$id),]
ancestry_index <- match(rownames(adult$meta),adult$ancestry$id)
adult$meta$ancestry <- adult$ancestry$ancestry[ancestry_index]

##Sort metadata by age, convert to days post conception log scale
#fetal convert to log pcd
colnames(fetal$meta) <- c("Subject", "age", "sex", "inferSex", "trimester", "ancestry", "study", "pcw")
fetal$meta$logPcd <- log(fetal$meta$pcw*7) #log of post conception days
fetal$meta <- fetal$meta[order(fetal$meta$logPcd), ]
fetal$eur_meta <- subset(fetal$meta, ancestry == "eur")
#get adult age in log pcd
adult$meta$age <- as.numeric(adult$meta$age)
adult$meta$logPcd <- log(adult$meta$age*365) #log of post conception day
adult$meta <- adult$meta[order(adult$meta$logPcd), ]
adult$meta$Subject <- rownames(adult$meta)
adult$eur_meta <- subset(adult$meta, ancestry == "EUR")
adult$eur_meta$ancestry <- tolower(adult$eur_meta$ancestry)
#combine metadata
#change ids to match the merge ids from genExp merge
fetal$eur_meta$Subject <- fetal$eur_meta$Subject %>% gsub("^(\\d)", "X\\1", .)
adult$eur_meta$Subject[adult$eur_meta$Subject %in% fetal$eur_meta$Subject] <- paste0(adult$eur_meta$Subject, ".a")[adult$eur_meta$Subject %in% fetal$eur_meta$Subject]
fetal$eur_meta$Subject[fetal$eur_meta$Subject %in% adult$eur_meta$Subject] <- paste0(fetal$eur_meta$Subject, ".f")[fetal$eur_meta$Subject %in% adult$eur_meta$Subject]
combined_meta <- merge(fetal$eur_meta, adult$eur_meta, all = T)
write.table(combined_meta, file = paste0(output_dir,"combo_meta.tsv"), row.names = FALSE, sep = "\t")