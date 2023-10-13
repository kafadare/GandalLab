##This is the code used to combine adult and fetal metadata, age is NOT converted to shared scale Sep 2023 EK
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


#reformat adult metadata for batch
adult$raw_meta <- as.data.frame(t(adult$raw_meta))
colnames(adult$raw_meta) <- adult$raw_meta[1, ]
adult$raw_meta <- adult$raw_meta[-1, ]
#add adult data batch column in metadata df
batch_index <- match(colData(adult_se)$names,rownames(adult$raw_meta))
adult$raw_meta$study <- colData(adult_se)$batch[batch_index]

#get index match for ancestry data, append column to meta data
adult_ancestry$id <- id_format_fix(adult_ancestry$id)
missing_ids <- rownames(adult$raw_meta)[!(rownames(adult$raw_meta) %in% adult_ancestry$id)] # 61 ids "missing" from the ancestry file
adult$raw_meta <- adult$raw_meta[(rownames(adult$raw_meta) %in% adult_ancestry$id),]
ancestry_index <- match(rownames(adult$raw_meta),adult_ancestry$id)
adult$raw_meta$ancestry <- adult_ancestry$ancestry[ancestry_index]

#Sort metadata by age
colnames(fetal$raw_meta) <- c("Subject", "age", "sex", "inferSex", "trimester", "ancestry", "study", "pcw")
fetal$raw_meta <- fetal$raw_meta[order(fetal$raw_meta$age), ]
fetal$eur_meta <- subset(fetal$raw_meta, ancestry == "eur")
adult$raw_meta$age <- as.numeric(adult$raw_meta$age)
adult$raw_meta <- adult$raw_meta[order(adult$raw_meta$age), ]
adult$raw_meta$Subject <- rownames(adult$raw_meta)
adult$eur_meta <- subset(adult$raw_meta, ancestry == "EUR")
adult$eur_meta$ancestry <- tolower(adult$eur_meta$ancestry)
#combine metadata
#change ids to match the merge ids from genExp merge
fetal$eur_meta$Subject <- fetal$eur_meta$Subject %>% gsub("^(\\d)", "X\\1", .)
fetal$eur_meta$Subject[fetal$eur_meta$Subject %in% adult$eur_meta$Subject] <- paste0(fetal$eur_meta$Subject, ".f")[fetal$eur_meta$Subject %in% adult$eur_meta$Subject]
adult$eur_meta$Subject[adult$eur_meta$Subject %in% fetal$eur_meta$Subject] <- paste0(adult$eur_meta$Subject, ".a")[adult$eur_meta$Subject %in% fetal$eur_meta$Subject]
combined_meta <- merge(fetal$eur_meta, adult$eur_meta, all = T)
write.table(combined_meta, file = paste0(output_dir,"combo_meta.tsv"), row.names = FALSE, sep = "\t")