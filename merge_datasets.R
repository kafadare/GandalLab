##Combine geneExp datasets and metadata for fetal and adult EUR

#Set up packages and functions
setwd("~/project-gandalm/GandalLab")
output_dir <- ("~/project-gandalm/comb_data/")
source("data_functions.R")
source("analysis_fcts.R")
# List of required packages
required_packages <- c("magrittr", "genio", "dplyr", "BiocManager", "snpStats", "SummarizedExperiment")
# Install or load missing packages
load_install_pkg(required_packages)

adult_se <- readRDS("/u/project/gandalm/kafadare/AdultBigBrain_gene_exp_raw_042923.RDS")

fetal$genExp.counts <- read.table("/u/project/gandalm/cindywen/isoform_twas/salmon/expression.final/gene.noVersion.scaled.counts.tsv", header = TRUE, sep = "\t") #counts file
fetal$genExp.tpm <- read.table("/u/project/gandalm/cindywen/isoform_twas/salmon/expression.final/gene.noVersion.TPM.tsv", header = TRUE, sep = "\t") #tpm file
fetal$raw_meta <- "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/metadata_654.tsv"

adult_se <- readRDS("/u/project/gandalm/kafadare/AdultBigBrain_gene_exp_raw_042923.RDS")
assays(adult_se) 
adult$raw_data.counts <- adult_se %>% assay(adult_struct)[1]
adult$raw_data.tpm <-  adult_se %>% assay(adult_struct)[2]
write.table(adult.counts,file="adult.counts.scaled.tsv",quote=FALSE, sep='\t')
write.table(adult.tpm,file="adult.TPM.tsv",quote=FALSE, sep='\t')
adult$raw_meta <- read.table("/u/project/gandalm/kafadare/cov_hcp0_gene.txt", header = TRUE, sep = "\t")
adult_ancestry <-  read.table("/u/project/gandalm/kafadare/pops.txt", header = F, sep = "\t", col.names = c("id", "ancestry"))

#fetal <- load_genExp_data(fetal_path, fetal_meta_path)
#adult <- load_genExp_data(adult_path, adult_meta_path)

#assay adult raw data gene counts
#adult$assay <- assay(adult$raw_data)
#reformat adult metadata
adult$raw_meta <- as.data.frame(t(adult$raw_meta))
colnames(adult$raw_meta) <- adult$raw_meta[1, ]
adult$raw_meta <- adult$raw_meta[-1, ]

#add adult data batch column in metadata df
colData(adult$raw_data)$names <- id_format_fix(colData(adult$raw_data)$names)
batch_index <- match(colData(adult$raw_data)$names,rownames(adult$raw_meta))
adult$raw_meta$study <- colData(adult$raw_data)$batch[batch_index]

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
#write.table(combined_meta, file = paste0(output_dir,"combo_meta.tsv"), row.names = FALSE, sep = "\t")

#Sort both genExp data by age
#clean up fetal data
#Add "X" in front of the ids that start with a number to match with the naming scheme in the geneExp file
f_ids_by_age <- fetal$eur_meta$Subject %>% intersect(.,colnames(fetal$raw_data))
fetal$sorted_data <- fetal$raw_data[,f_ids_by_age]
#sort adult data
colnames(adult$assay) <- id_format_fix(colnames(adult$assay))
a_ids_by_age <- rownames(adult$eur_meta) %>% intersect(.,colnames(adult$assay))
adult$sorted_data <- adult$assay[,a_ids_by_age]

#remove gencode version number from the adult dataset gene id
rownames(adult$sorted_data) <- substring(rownames(adult$sorted_data), 1, 15)

#get vector of shared genes and subset both fetal and adult data for NON-shared genes
shared_genes <- rownames(adult$sorted_data)[which(rownames(adult$sorted_data) %in% rownames(fetal$sorted_data))]
fetal_only_genes <- fetal$sorted_data[which(!(rownames(fetal$sorted_data) %in% rownames(adult$sorted_data))),]
fetal_only_genes$genid <- rownames(fetal_only_genes)
adult_only_genes <- as.data.frame(adult$sorted_data[which(!(rownames(adult$sorted_data) %in% rownames(fetal$sorted_data))),])
adult_only_genes$genid <- rownames(adult_only_genes)

#combine the two datasets for shared genes
combined_genExp <- merge(fetal$sorted_data, adult$sorted_data, by = "row.names", all = F, suffixes = c(".f", ".a"))
colnames(combined_genExp)[colnames(combined_genExp) == 'Row.names'] <- 'genid'
#write.csv(my_dataframe, file = "output_file.csv", row.names = FALSE)

#get Chr number and start and end codon
bm_gene_fetal <- get_gene_info(fetal_only_genes$genid)
bm_gene_adult <- get_gene_info(adult_only_genes$genid)
bm_gene_combo <- get_gene_info(combined_genExp$genid)

#add gene info to the dataset
combined_genExp <- merge(bm_gene_combo, combined_genExp, by = "genid")
fetal_only_genes <- merge(bm_gene_fetal, fetal_only_genes, by = "genid")
adult_only_genes <- merge(bm_gene_adult, adult_only_genes, by = "genid")

##should save again!!!
write.table(combined_genExp, file = paste0(output_dir,"combo_genExp.tsv"), row.names = FALSE, sep = "\t")
write.table(fetal_only_genes, file = paste0(output_dir,"fetal_only_genExp.tsv"), row.names = FALSE, sep = "\t")
write.table(adult_only_genes, file = paste0(output_dir,"adult_only_genExp.tsv"), row.names = FALSE, sep = "\t")

#coldata listdata batch for batch, bampaths are also there
#fetal_init_cols <- fetal$raw_data[,2:4] #removed first column because first and second column are identical
#fetal_subj_cols <- fetal$raw_data[,5:ncol(fetal$raw_data)]