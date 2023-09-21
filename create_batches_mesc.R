#vcf_to_plink <- function (vcf_path, geneExp_path, cov_file_path, batch = FALSE, batch_size, batch_step) {
#load functions from data_functions.R
setwd("~/project-gandalm/GandalLab")
output_dir <- ("/u/scratch")
source("data_functions.R")
# List of required packages
required_packages <- c("magrittr", "genio", "dplyr", "BiocManager", "snpStats")
# Install or load missing packages
load_install_pkg(required_packages)
#Load in data
#bed_path <- "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeRel.bed"
#bim_path <- "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeRel.bim"
#fam_path <- "/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeRel.fam"
#plink_data <- snpStats::read.plink(bed_path, bim_path, fam_path)
#bim_data <- read_plink(bim_path)
#fam_data <- read_plink(fam_path)

geneExp_path <- "/u/project/gandalm/cindywen/isoform_twas/TWAS/data/eur_gene_exp_regressed.txt"
geneExp <- read.table(geneExp_path, header = TRUE, sep = "\t")

metadata_path <- "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/metadata_654.tsv"
metadata <- read.table(metadata_path, header = TRUE, sep = "\t") %>% filter(ancestry == "eur") %>% arrange(Age)
##Question: metadata filtered for EUR has 292 samples, vcf file has 281 samples. GeneExp has 284 samples. Where does the discrepancy come from?

#Match gene exp and gwas data with metadata age
#sorted_genotype <- plink_data$genotypes[ids_by_age,]
#sorted_fam <- plink_data$fam(match)
# Create a new data frame with matched rows
geneExp_init_cols <- geneExp[,1:4]
geneExp_subj_cols <- geneExp[,5:ncol(geneExp)]
#Add "X" in front of the ids that start with a number to match with the naming scheme in the geneExp file
ids_by_age <- metadata$Subject %>% gsub("^(\\d)", "X\\1", .) %>% intersect(.,colnames(geneExp_subj_cols))
geneExp_subj_cols_sorted <- geneExp_subj_cols[,ids_by_age]

##Create batches by column
#Define batchsize and batchstep
batchsize = 150
batchstep = 25

geneExp_sorted <- cbind(geneExp_init_cols, geneExp_subj_cols_sorted)

  
  #Create the batches
  create_batch()
}
#save batches as vcf in scratch