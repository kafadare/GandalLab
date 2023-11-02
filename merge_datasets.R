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

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
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
fetal$eur_meta$Subject <- id_format_fix(fetal$eur_meta$Subject)
adult$eur_meta$Subject[adult$eur_meta$Subject %in% fetal$eur_meta$Subject] <- paste0(adult$eur_meta$Subject, ".a")[adult$eur_meta$Subject %in% fetal$eur_meta$Subject]
combined_meta <- merge(fetal$eur_meta, adult$eur_meta, all = T)
write.table(combined_meta, file = paste0(output_dir,"combo_meta.tsv"), row.names = FALSE, sep = "\t")

#Sort both genExp data by age, both counts and tpm
f_ids_by_age <- fetal$eur_meta$Subject %>% intersect(.,colnames(fetal$counts))
fetal$sorted.counts<- fetal$counts[,f_ids_by_age]
f_ids_by_age <- fetal$eur_meta$Subject %>% intersect(.,colnames(fetal$tpm))
fetal$sorted.tpm <- fetal$tpm[,f_ids_by_age]
#write.table(fetal$sorted.counts, file = paste0(output_dir,"fetal_eur_counts.tsv"), row.names = FALSE, sep = "\t")
#write.table(fetal$sorted.tpm, file = paste0(output_dir,"fetal_eur_tpm.tsv"), row.names = FALSE, sep = "\t")

#sort adult data
colnames(adult$counts) <- id_format_fix(colnames(adult$counts))
a_ids_by_age <- rownames(adult$eur_meta) %>% intersect(.,colnames(adult$counts))
adult$sorted.counts <- adult$counts[,a_ids_by_age]
colnames(adult$tpm) <- id_format_fix(colnames(adult$tpm))
a_ids_by_age <- rownames(adult$eur_meta) %>% intersect(.,colnames(adult$tpm))
adult$sorted.tpm <- adult$tpm[,a_ids_by_age]

#remove gencode version number from the adult dataset gene id
rownames(adult$sorted.counts) <- substring(rownames(adult$sorted.counts), 1, 15)
rownames(adult$sorted.tpm) <- substring(rownames(adult$sorted.tpm), 1, 15)

#get vector of shared genes and subset both fetal and adult data for NON-shared genes
shared_genes <- rownames(adult$sorted.counts)[which(rownames(adult$sorted.counts) %in% rownames(fetal$sorted.counts))]
fetal_only_genes <- fetal$sorted.counts[which(!(rownames(fetal$sorted.counts) %in% rownames(adult$sorted.counts))),]
fetal_only_genes$genid <- rownames(fetal_only_genes)
adult_only_genes <- as.data.frame(adult$sorted.counts[which(!(rownames(adult$sorted.counts) %in% rownames(fetal$sorted.counts))),])
adult_only_genes$genid <- rownames(adult_only_genes)

#combine the two datasets for shared genes
comb_genExp.counts <- merge(fetal$sorted.counts, adult$sorted.counts, by = "row.names", all = F, suffixes = c(".f", ".a"))
colnames(comb_genExp.counts)[colnames(comb_genExp.counts) == 'Row.names'] <- 'genid'
comb_genExp.tpm <- merge(fetal$sorted.tpm, adult$sorted.tpm, by = "row.names", all = F, suffixes = c(".f", ".a"))
colnames(comb_genExp.tpm)[colnames(comb_genExp.tpm) == 'Row.names'] <- 'genid'
#write.csv(my_dataframe, file = "output_file.csv", row.names = FALSE)'

#get Chr number and start and end codon
gtf_fetal <- get_gene_info(fetal_only_genes$genid)
gtf_adult <- get_gene_info(adult_only_genes$genid)
gtf_comb <- get_gene_info(comb_genExp.counts$genid)

#save these gtf files too!! ^^^
write.table(gtf_comb, file = paste0(output_dir,"gene_info_comb.tsv"), row.names = FALSE, sep = "\t")
write.table(gtf_fetal, file = paste0(output_dir,"gene_info_fetal.tsv"), row.names = FALSE, sep = "\t")
write.table(gtf_adult, file = paste0(output_dir,"gene_info_adult.tsv"), row.names = FALSE, sep = "\t")

#add gene info to the dataset
#comb_genExp.counts <- merge(gtf_comb, combined_genExp, by = "genid")
#comb_genExp.counts <- merge(gtf_comb, combined_genExp, by = "genid")
#fetal_only_genes <- merge(gtf_fetal, fetal_only_genes, by = "genid")
#adult_only_genes <- merge(gtf_adult, adult_only_genes, by = "genid")

##should save again!!!
write.table(comb_genExp.counts, file = paste0(output_dir,"combo_genExp_counts.tsv"), row.names = FALSE, sep = "\t")
write.table(comb_genExp.tpm, file = paste0(output_dir,"combo_genExp_tpm.tsv"), row.names = FALSE, sep = "\t")
write.table(fetal_only_genes, file = paste0(output_dir,"fetal_only_genExp.tsv"), row.names = FALSE, sep = "\t")
write.table(adult_only_genes, file = paste0(output_dir,"adult_only_genExp.tsv"), row.names = FALSE, sep = "\t")

#coldata listdata batch for batch, bampaths are also there
#fetal_init_cols <- fetal$raw_data[,2:4] #removed first column because first and second column are identical
#fetal_subj_cols <- fetal$raw_data[,5:ncol(fetal$raw_data)]
