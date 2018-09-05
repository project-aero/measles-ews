#!/bin/bash
#SBATCH -t 20:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 20
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
 
cd $SBATCH_O_WORKDIR
module load R
 
# parallel R: submit job with one MPI master
R CMD BATCH global-search-mif.R