#!/bin/bash -l
#$ -cwd
#$ -l h_data=4G,h_rt=24:00:00
#$ -j y
#$ -o /u/home/k/kafadare/project-gandalm/mesc_genExp_scores/
#$ -m a
#$ -t 1-2

CHR=${SGE_TASK_ID}

geno="/u/home/k/kafadare/project-gandalm/merged_genotype_Chr/merged_chr${CHR}" ##file determined by which chromosome

keep_file="/u/home/k/kafadare/project-gandalm/comb_data/new_data/processed/genoID_1968.txt"

outdir="/u/home/k/kafadare/project-gandalm/merged_genotype_Chr/merged_chr_filtered${CHR}"

. /u/local/Modules/default/init/modules.sh
module load plink

plink --bfile ${geno} --keep ${keep_file} --make-bed --out ${outdir}


