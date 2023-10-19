setwd("~/project-gandalm/GandalLab")
output_dir <- ("~/project-gandalm/comb_data/processed")
plot_dir <- "/u/home/k/kafadare/project-gandalm/plots"
source("data_functions.R")
source("analysis_fcts.R")
# List of required packages
required_packages <- c("magrittr", "dplyr", "data.table", "DESeq2", "WGCNA", "sva", "ggplot2")
#required_packages <- c("magrittr", "dplyr", "data.table", "ggplot2")
# Install or load missing packages
load_install_pkg(required_packages)
#library(DESeq2)
#library(sva)
#library(WGCNA)

#Load Data
datExpr <- fread("/u/home/k/kafadare/project-gandalm/comb_data/new_data/combo_genExp_counts.tsv", header = T, data.table = F)
datExpr.tpm <- fread("/u/home/k/kafadare/project-gandalm/comb_data/new_data/combo_genExp_tpm.tsv", header = T, data.table = F)
gtf <- fread("/u/home/k/kafadare/project-gandalm/comb_data/new_data/gene_info_comb.tsv", header = T, data.table = F)
#colnames(gtf_comb)[2] <- "Chr"
meta <- fread("/u/home/k/kafadare/project-gandalm/comb_data/combo_meta.tsv", header = T, data.table = F)
##Change names in meta file to match the names in the genExp files (especially for the overlap samples)

# test<- prep_datExp(datExpr, gtf) 
overlap_ids <- meta[grep(".f",meta$Subject),'Subject']
overlap_ids <- gsub(".f", "", overlap_ids)
ind <- which(meta$Subject %in% overlap_ids)
meta[ind, 1] <- gsub("(\\d)$", "\\1.a", meta[ind, 1])
# 
# data.batch <- c()
# for (i in 1:ncol(test)) {
#   sample <- colnames(test)[i]
#   #match the sample to id in meta file and find corresponding "study"
#   data.batch[i] <- meta[which(meta$Subject %in% sample), "study"]
# }


norm_batch <- prep_datExp(datExpr, gtf) %>% norm_batch(., meta)
datExpr <- norm_batch$datExpr
combat_expr <- norm_batch$combat

# Save the plot
pdf(paste0(plot_dir,"ComBat_processed_shared.pdf"))
plot_norm <- prep_datExp(combat_expr, gtf_comb) %>% plot_gen_density(., title_str = "Norm & Batch Correct", ylim = c(0, 1))
dev.off() 

# phenotype.bed
gtf <- gtf %>% filter(V5 %in% rownames(datExpr.final))

# prepare BED file
setDT(combat_expr, keep.rownames = TRUE)
bed <- merge(combat_expr, gtf, by.x = "rn", by.y = "V5", all = TRUE)
bed <- as.data.frame(bed)

bed <- bed[, c(ncol(bed)-3, ncol(bed)-2, ncol(bed)-1, ncol(bed), 1:(ncol(bed)-4))]

# special note! strand; 0/1-based conversion
bed$V1 <- gsub("chr","",bed$V1)

for (i in 1:dim(bed)[1]){
  if (bed[i,4] == "+") {
    gtf_start=bed[i,2]
    bed[i,2]=gtf_start-1
    bed[i,3]=gtf_start
  }
  else if (bed[i,4] == "-") {
    gtf_end=bed[i,3]
    bed[i,2]=gtf_end-1
    bed[i,3]=gtf_end
  }
}

bed <- bed[order(bed$V1, bed$V2),]

bed <- bed %>% select(-V4)
colnames(bed)[1:4] <- c("Chr", "start", "end", "ID")

write.table(bed, paste0(output_dir, ".counts.scaled.normalized.bed"), quote=F, sep="\t", row.names=F, col.names=T)


