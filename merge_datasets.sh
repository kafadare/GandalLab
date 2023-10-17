### R_submit.sh START ###
#!/bin/bash
#$ -cwd
#$ -o  joblog.$JOB_ID
#$ -j y
#  Resources requested
#  PLEASE CHANGE THE RESOURCES REQUESTED AS NEEDED:
#$ -l h_data=4G,h_rt=1:00:00
#  PLEASE CHANGE THE NUMBER OF CORES REQUESTED AS NEEDED:
#$ -pe shared 1
#  Email address to notify
#$ -M $USER@mail
#$ -m bea
# #$ -V

#
# Output job info on joblog file:
#
echo " "
echo "Job myscript, ID no. $JOB_ID started on:   "` hostname -s `
echo "Job myscript, ID no. $JOB_ID started at:   "` date `
echo " "

#
# Set up job environment:
#
. /u/local/Modules/default/init/modules.sh
module load R
# SHOULD YOU NEED TO LINK TO ANY PARTICULAR LIBRARY, USE:
# export LD_LIBRARY_PATH=/path/to/libdir/if/needed:$LD_LIBRARY_PATH
# SHOULD YOU NEED TO CHANGE THE LOCATION OF THE R TEMPORARY DIRECTORY:
# if [ ! -d $SCRATCH/R_scratch ]; then
#    mkdir $SCRATCH/R_scratch
# fi
# export TMPDIR=$SCRATCH

#
# Run the R script:
#
echo " "
#  (SUBSTITUTE THE NAME OF YOUR R SCRIPT AND OUTPUT BELOW):
echo R CMD BATCH --no-save --no-restore  myscript.R output.$JOB_ID
#  (SUBSTITUTE THE NAME OF YOUR R SCRIPT AND OUTPUT BELOW):
/usr/bin/time -v R CMD BATCH --no-save --no-restore  myscript.R  output.$JOB_ID
#
echo " "
echo "Job myscript, ID no. $JOB_ID finished at:  "` date `
echo " "
### R_submt.sh STOP ###