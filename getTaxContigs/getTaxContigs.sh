#!/bin/sh

#SBATCH --job-name="getTaxa"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH -c 1

### adjust/add sbatch flags as needed

cd $SLURM_SUBMIT_DIR

module load Biopython/1.78-foss-2020b 
mkdir refs_for_annotation

python /home/aknyshov/alex_data/tree_shew_analysis/TreeshrewProject/getTaxContigs/getTaxContigs.py reftaxa.txt aligned_loci refs_for_annotation