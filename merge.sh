#!/bin/bash -lj
#$ -cwd
#$ -j y
#$ -sync y

module load bcftools
bcftools isec -p isec adult_EUR.vcf.gz /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeGeneOutlier.vcf.gz
bcftools merge -o merged.vcf adult_EUR.vcf.gz /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeGeneOutlier.vcf.gz
