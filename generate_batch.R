##Script to generate batches for sliding window analysis

setwd("/u/project/gandalm/kafadare/GandalLab")
batch_output_dir <- ("/u/project/gandalm/kafadare/comb_data/batches/")
plot_dir <- ("/u/project/gandalm/kafadare/plots/")
source("data_functions.R")
source("analysis_fcts.R")
# List of required packages
required_packages <- c("magrittr", "dplyr", "data.table", "ggplot2")
# Install or load missing packages
load_install_pkg(required_packages)

#change the datExpr file location after regressing out covs and/or removing the data without genotypes & matching IDs with genotypes
datExpr <- fread("/u/home/k/kafadare/project-gandalm/comb_data/new_data/processed/genExp_covRegOut.tsv", data.table=F)
datExpr <- datExpr %>% dplyr::rename(., genid = V1)
gtf <-  fread("/u/home/k/kafadare/project-gandalm/comb_data/new_data/gene_info_comb.tsv", data.table=F)
colnames(gtf)
#change covs file location after saving full covs file with PCs and genoPCs
covs <- fread("/u/home/k/kafadare/project-gandalm/comb_data/new_data/combined_meta.tsv", header = T, data.table = F)

#MESC format 
# Col 1 GENE: Gene name
# Col 2 CHR: Chromosome number of gene
# Col 3 GENE_COORD: Chromosomal coordinate of gene start or midpoint. Because we take a 1MB region around this point, the exact location isn't too important.
# 4th column onwards: Gene expression values with column name corresponding to sample name


# prepare BED file
setDT(datExpr, keep.rownames = F)
mescGE <- merge(datExpr, gtf, by = "genid", all = TRUE)
mescGE <- as.data.frame(mescGE)

mescGE <- mescGE[, c(1, ncol(mescGE)-3, ncol(mescGE)-2, ncol(mescGE)-1, ncol(mescGE), 2:(ncol(mescGE)-4))]
head(colnames(mescGE))

# special note! strand; 0/1-based conversion
for (i in 1:dim(mescGE)[1]){
  if (mescGE[i,4] == "+") {
    gtf_start=mescGE[i,2]
    mescGE[i,2]=gtf_start-1
    mescGE[i,3]=gtf_start
  }
  else if (mescGE[i,4] == "-") {
    gtf_end=mescGE[i,3]
    mescGE[i,2]=gtf_end-1
    mescGE[i,3]=gtf_end
  }
}

mescGE <- mescGE[order(mescGE$Chr, mescGE$start),]
#remove MT from genExp matrix
mescGE <- mescGE %>% subset(Chr != "MT")
mescGE$strand <- NULL
mescGE$end <- NULL
head(colnames(mescGE))
colnames(mescGE)[1:3] <- c("GENE", "CHR", "GENE_COORD")
head(colnames(mescGE))

batches <- create_batch(mescGE, batch_size = 150, batch_step = 25, method = "sliding_window") %>% save_list_elements(., delimiter = "/t", file_extension = "tsv", save_folder = batch_output_dir, file_name = "GE_batch")
#last batch sample n = 133
#proof of concept, check, but need to figure out what is the actual data we are using.


