setwd("~/project-gandalm/GandalLab")
output_dir <- ("~/project-gandalm/comb_data/processed/")
plot_dir <- "/u/home/k/kafadare/project-gandalm/plots/"
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

datExpr <- prep_datExp(datExpr, gtf)

norm_batch <- prep_datExp(datExpr, gtf) %>% norm_batch(., meta, output_dir)
datExpr <- norm_batch$datExpr
combat_expr <- norm_batch$combat

datExpr.final <- fread("/u/home/k/kafadare/project-gandalm/comb_data/processed/counts.norm.noComBat.tsv", data.table = F)
rownames(datExpr.final) <- datExpr.final$V1
datExpr.final <- datExpr.final[,-1]
combat_expr <- fread("/u/home/k/kafadare/project-gandalm/comb_data/processed/counts.batch.processed.tsv", data.table = F)
rownames(combat_expr) <- combat_expr$V1
combat_expr <- combat_expr[,-1]

# Normalize
datExpr.vst <- varianceStabilizingTransformation(as.matrix(datExpr), blind = TRUE)
datExpr.vst <- as.data.frame(datExpr.vst)
#### 4
# Remove outliers
normadj <- adjacency(datExpr.vst,type = 'signed',corFnc = 'bicor')   #Calculate network adjacency
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity   #Extract connectivity of each sample
Z.C <- (C-mean(C))/sqrt(var(C))   #Covert to Z-score
outliers <- (Z.C < -3)
outlier.df <- data.frame(c(which(outliers)))
outlier.df$c.which.outliers.. <- rownames(outlier.df)
write.table(outlier.df,paste0(output_dir, ".outlier.txt"), sep="\t", quote=F, col.names = F, row.names = F)

# phenotype.bed
gtf <- gtf %>% filter(genid %in% rownames(datExpr.final))

# prepare BED file
setDT(combat_expr, keep.rownames = TRUE)
bed <- merge(combat_expr, gtf, by.x = "rn", by.y = "genid", all = TRUE)
bed <- as.data.frame(bed)

bed <- bed[, c( ncol(bed) -3, ncol(bed)-2, ncol(bed)-1, ncol(bed), 1:(ncol(bed)-3))]

# special note! strand; 0/1-based conversion
bed$Chr <- gsub("chr","",bed$Chr)

#see what the +/- is
for (i in 1:dim(bed)[1]){
  if (bed[i,4] == +1) {
    gtf_start=bed[i,2]
    bed[i,2]=gtf_start-1
    bed[i,3]=gtf_start
  }
  else if (bed[i,4] == -1) {
    gtf_end=bed[i,3]
    bed[i,2]=gtf_end-1
    bed[i,3]=gtf_end
  }
}

bed <- bed[order(bed$Chr, bed$start),]

bed <- bed %>% select(-strand)
colnames(bed)[1:4] <- c("Chr", "start", "end", "ID")

write.table(bed, paste0(output_dir, "counts.scaled.normalized.bed"), quote=F, sep="\t", row.names=F, col.names=T)


