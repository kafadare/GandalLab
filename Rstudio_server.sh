#!/bin/bash
space="$1" #in GB
time="$2" #in hours
qrsh -l h_data=$dataG,h_rt=$time:00:00
# Create small tmp directories for RStudio to write into
mkdir -pv $SCRATCH/rstudiotmp/var/lib
mkdir -pv $SCRATCH/rstudiotmp/var/run
mkdir -pv $SCRATCH/rstudiotmp/tmp
#Setup apptainer
module load apptainer
#Run rstudio
apptainer run -B $SCRATCH/rstudiotmp/var/lib:/var/lib/rstudio-server -B $SCRATCH/rstudiotmp/var/run:/var/run/rstudio-server -B $SCRATCH/rstudiotmp/tmp:/tmp $H2_CONTAINER_LOC/h2-rstudio_4.1.0.sif
# This command will display some information and a `ssh -L ...` command for you to run on a separate terminal
