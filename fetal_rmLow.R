##Code to remove lowly expressed genes for fetal genExp data, and save plots for before and after removal

#Set up packages and functions
setwd("~/project-gandalm/GandalLab")
output_dir <- ("~/project-gandalm/")
plot_dir <- ("~/project-gandalm/plots/")
source("data_functions.R")
source("analysis_fcts.R")
# List of required packages
#required_packages <- c("magrittr", "genio", "dplyr", "BiocManager", "SummarizedExperiment")
required_packages <- c("magrittr","dplyr")
# Install or load missing packages
#load_install_pkg(required_packages)
library(magrittr)
library(dplyr)

fetal_counts_path <- "/u/project/gandalm/cindywen/isoform_twas/salmon/expression.final/gene.noVersion.scaled.counts.tsv"
fetal_tpm_path <- "/u/project/gandalm/cindywen/isoform_twas/salmon/expression.final/gene.noVersion.TPM.tsv"
fetal <- load_data(names = c("counts", "tpm"), fetal_counts_path, fetal_tpm_path)
fetal$counts$genid <- rownames(fetal$counts)
fetal$tpm$genid <- rownames(fetal$tpm)
gtf_fetal <- get_gene_info(fetal$counts$genid)
#All fetal genes
# Save the plot
pdf(paste0(plot_dir,"fetal_all_counts.pdf"))
# 2. Create a plot
all_plot_counts <- prep_datExp(fetal$counts, gtf_fetal) %>% plot_gen_density(., title_str = "Fetal Counts All", ylim = c(0, 0.75))
# Close the pdf file
dev.off() 
# Save the plot
pdf(paste0(plot_dir,"fetal_all_tpm.pdf"))
# 2. Create a plot
all_plot_tpm <- prep_datExp(fetal$tpm, gtf_fetal) %>% plot_gen_density(., title_str = "Fetal Tpm All", ylim = c(0, 2.0))
# Close the pdf file
dev.off() 

#Remove lowly expressed genes <0.1 25%
fetal_lowrm <- prep_datExp(fetal$counts, gtf_fetal) %>% rm_low(., fetal$tpm, gtf_fetal, cutoff = 0.1, percent = 0.25)
fetal_lowrm_tpm <- prep_datExp(fetal$tpm, gtf_fetal) %>% rm_low(., fetal$tpm, gtf_fetal, cutoff = 0.1, percent = 0.25)

print(cat("Number of genes before removal:", dim(fetal$counts), "Number of genes after removing lowly expressed:",dim(fetal_lowrm), sep='\n'))

#save the data with lowly expressed removed
write.table(fetal_lowrm,file = paste0(output_dir,"fetal_rmlow_counts.tsv"),quote=FALSE, sep='\t')
write.table(fetal_lowrm_tpm,file = paste0(output_dir,"fetal_rmlow_tpm.tsv"),quote=FALSE, sep='\t')

# Save the plot after removal
pdf(paste0(plot_dir,"fetal_rmLow_counts.pdf"))
# 2. Create a plot
rmLow_plot_counts <-  plot_gen_density(fetal_lowrm, title_str = "Fetal Low Exp Rm", ylim = c(0, 0.75))
# Close the pdf file
dev.off() 


