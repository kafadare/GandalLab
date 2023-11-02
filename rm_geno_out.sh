#!/bin/bash -l
#$ -cwd
#$ -l h_data=XXX,h_rt=24:00:00
#$ -j y
#$ -o /u/home/k/kafadare/project-gandalm/mesc_genExp_scores/
#$ -m a
#$ -t 1-22

CHR=${SGE_TASK_ID}

expr = /u/home/k/kafadare/project-gandalm/genExp_batches

geno=/u/home/k/kafadare/project-gandalm/merged_genotype_Chr/XXX ##file determined by which chromosome

cov=/u/home/k/kafadare/project-gandalm/comb_data/new_data/processed/

rel_file=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/related.txt

outdir=/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_${num_hcp}hcp

