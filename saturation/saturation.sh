#!/bin/bash
#SBATCH --job-name="satur"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G

aligned_loci_path="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/aligned_loci/"
script_path="/home/aknyshov/alex_data/tree_shew_analysis/TreeshrewProject/saturation/saturation.R"
cd $SLURM_SUBMIT_DIR
module load R/4.0.3-foss-2020b
date
mkdir saturation_assessments
Rscript ${script_path} ${aligned_loci_path}
mv saturation_output.csv saturation_assessments/
date