#!/bin/bash

#SBATCH --job-name=CePTER
#SBATCH --partition=fuchs
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=1000   
#SBATCH --no-requeue 
#SBATCH --mail-type=ALL

R < 00_Installer.R --no-save --slave
R < 01_Mapping_Alignment.R --no-save --slave
