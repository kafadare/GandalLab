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

adult_se <- readRDS("/u/project/gandalm/kafadare/AdultBigBrain_tx_exp_raw_042923.RDS")
#gene_raw = tximport::summarizeToGene(adult_se)
fetal <- list()
fetal$genExp.counts <- read.table("/u/project/gandalm/cindywen/isoform_twas/salmon/expression.final/gene.noVersion.scaled.counts.tsv", header = TRUE, sep = "\t") #counts file
fetal$genExp.tpm <- read.table("/u/project/gandalm/cindywen/isoform_twas/salmon/expression.final/gene.noVersion.TPM.tsv", header = TRUE, sep = "\t") #tpm file
fetal$raw_meta <- read.table("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/metadata_654.tsv", header = TRUE, sep = "\t") 

adult_se <- readRDS("/u/home/k/kafadare/project-gandalm/AdultBigBrain_gene_exp_raw_042923.RDS")
assays(adult_se)
adult <- list()
adult$genExp.counts <- adult_se %>% assay(.,1)
adult$genExp.tpm <-  adult_se %>% assay(.,2)
#write.table(adult$raw_data.counts,file="adult.counts.scaled.tsv",quote=FALSE, sep='\t')x
#write.table(adult$raw_data.tpm,file="adult.TPM.tsv",quote=FALSE, sep='\t')
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
colData(adult_se)$names <- id_format_fix(colData(adult_se)$names)
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
#write.table(combined_meta, file = paste0(output_dir,"combo_meta.tsv"), row.names = FALSE, sep = "\t")

#Sort both genExp data by age, both counts and tpm
f_ids_by_age <- fetal$eur_meta$Subject %>% intersect(.,colnames(fetal$genExp.counts))
fetal$sorted.counts<- fetal$genExp.counts[,f_ids_by_age]
f_ids_by_age <- fetal$eur_meta$Subject %>% intersect(.,colnames(fetal$genExp.tpm))
fetal$sorted.tpm <- fetal$genExp.tpm[,f_ids_by_age]
#write.table(fetal$sorted.counts, file = paste0(output_dir,"fetal_eur_counts.tsv"), row.names = FALSE, sep = "\t")
#write.table(fetal$sorted.tpm, file = paste0(output_dir,"fetal_eur_tpm.tsv"), row.names = FALSE, sep = "\t")

#sort adult data
colnames(adult$genExp.counts) <- id_format_fix(colnames(adult$genExp.counts))
a_ids_by_age <- rownames(adult$eur_meta) %>% intersect(.,colnames(adult$genExp.counts))
adult$sorted.counts <- adult$genExp.counts[,a_ids_by_age]
colnames(adult$genExp.tpm) <- id_format_fix(colnames(adult$genExp.tpm))
a_ids_by_age <- rownames(adult$eur_meta) %>% intersect(.,colnames(adult$genExp.tpm))
adult$sorted.tpm <- adult$genExp.tpm[,a_ids_by_age]

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
comb_genExp.counts <- merge(bm_gene_combo, combined_genExp, by = "genid")
comb_genExp.counts <- merge(bm_gene_combo, combined_genExp, by = "genid")
fetal_only_genes <- merge(bm_gene_fetal, fetal_only_genes, by = "genid")
adult_only_genes <- merge(bm_gene_adult, adult_only_genes, by = "genid")

##should save again!!!
write.table(comb_genExp.counts, file = paste0(output_dir,"combo_genExp_counts.tsv"), row.names = FALSE, sep = "\t")
write.table(comb_genExp.tpm, file = paste0(output_dir,"combo_genExp_tpm.tsv"), row.names = FALSE, sep = "\t")
write.table(fetal_only_genes, file = paste0(output_dir,"fetal_only_genExp.tsv"), row.names = FALSE, sep = "\t")
write.table(adult_only_genes, file = paste0(output_dir,"adult_only_genExp.tsv"), row.names = FALSE, sep = "\t")

#coldata listdata batch for batch, bampaths are also there
#fetal_init_cols <- fetal$raw_data[,2:4] #removed first column because first and second column are identical
#fetal_subj_cols <- fetal$raw_data[,5:ncol(fetal$raw_data)]