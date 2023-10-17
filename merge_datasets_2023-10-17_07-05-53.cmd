#!/bin/bash
#$ -cwd
#$ -o  merge_datasets.joblog.$JOB_ID
#$ -j y
#  Resources requested:
#$ -l h_data=32G,h_rt=2:00:00
#$ -pe shared 1
#  Email address to notify
#$ -M kafadare@mail
#$ -m bea
# #$ -V

#
# Output job info on joblog file:
#
echo " "
echo "Job merge_datasets, ID no. $JOB_ID started on:   "` hostname -s `
echo "Job merge_datasets, ID no. $JOB_ID started at:   "` date `
echo " "

#
# Set up job environment:
#
. /u/local/Modules/default/init/modules.sh
#module load R/4.2.2-BIO
 module load R/4.2.2-BIO
module list
#export LD_LIBRARY_PATH=/path/to/libdir/if/needed:$LD_LIBRARY_PATH   

# setting the number of threads to the number of cores requested:
export OMP_NUM_THREADS=1

#
# Run the R script:
#
echo ""
echo "R CMD BATCH --no-save --no-restore  merge_datasets.R merge_datasets.out.$JOB_ID"
#
/usr/bin/time -p R CMD BATCH --no-save --no-restore  merge_datasets.R  merge_datasets.out.$JOB_ID 
#
echo " "
echo "Job merge_datasets, ID no. $JOB_ID finished at:  "` date `
echo " "
