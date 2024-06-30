#!/bin/bash
#SBATCH -c 64
#SBATCH --partition=common,scavenger
#SBATCH --output="/hpc/group/hofflab/elb75/slurmout/R-%A_%a.out"
#SBATCH --error="/hpc/group/hofflab/elb75/slurmout/R-%A_%a.err"
#SBATCH --mem=228G # GB RAM
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=elb75@duke.edu
#SBATCH --exclude=dcc-ultrasound-04,dcc-ultrasound-06,dcc-ultrasound-07

module load R
## Export R_LIBS_USER = ~/R/x86_64-pc-linux-gnu-library/4.0
## R CMD BATCH SIMULATION.R
Rscript simulation.R
