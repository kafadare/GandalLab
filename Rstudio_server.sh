#!/bin/bash
echo $PATH
$PATH = "/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/bin/intel64:/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/bin:/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mpi/intel64/libfabric/bin:/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin:/u/local/compilers/intel/2020.4/debugger_2020/gdb/intel64/bin:/u/systems/UGE8.6.4/bin/lx-amd64:/u/local/bin:/u/local/sbin:/u/local/Modules/4.7.0/gcc-4.8.5/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/u/home/k/kafadare/bin"
space="$1" #in GB
time="$2" #in hours
qrsh -l h_data=${space},h_rt=${time}:00:00
# Create small tmp directories for RStudio to write into
mkdir -pv $SCRATCH/rstudiotmp/var/lib
mkdir -pv $SCRATCH/rstudiotmp/var/run
mkdir -pv $SCRATCH/rstudiotmp/tmp
#Setup apptainer
module load apptainer
#Run rstudio
apptainer run -B $SCRATCH/rstudiotmp/var/lib:/var/lib/rstudio-server -B $SCRATCH/rstudiotmp/var/run:/var/run/rstudio-server -B $SCRATCH/rstudiotmp/tmp:/tmp $H2_CONTAINER_LOC/h2-rstudio_4.1.0.sif
# This command will display some information and a `ssh -L ...` command for you to run on a separate terminal
