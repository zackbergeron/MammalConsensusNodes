#!/bin/bash
#SBATCH --job-name="screen"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G

concattree_path="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing_v2/sptrees.tre"
genetree_path="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing_v2/gtrees.tre"
genetree_names="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing_v2/gtrees.txt"
script_path="/home/aknyshov/alex_data/tree_shew_analysis/TreeshrewProject/screening/treescreen.R"

cd $SLURM_SUBMIT_DIR

module load R/4.0.3-foss-2020b

date
mkdir screening_assessments
Rscript ${script_path} ${concattree_path} ${genetree_path} ${genetree_names}
mv screening_output.csv screening_assessments/
date