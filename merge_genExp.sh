#!/bin/bash -l
#$ -cwd
#$ -j y
#$ -sync y

module load R
Rscript merge_datasets.R
