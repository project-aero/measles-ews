#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 21
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
 
cd ~/measles/code/
module load R
 
# parallel R: submit job with one MPI master
R CMD BATCH global-search-mif.R
