#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -sync y

module load R/4.3.0
Rscript merge_datasets.R
