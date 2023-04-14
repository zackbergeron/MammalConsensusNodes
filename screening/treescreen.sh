#!/bin/bash
#SBATCH --job-name="screen"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH --mail-user="zbergeron@uri.edu"
#SBATCH --mail-type=ALL
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G

concattree_path="/data/schwartzlab/zbergeron/TreeshrewProject/iqtree/sptrees.tre"
genetree_path="/data/schwartzlab/zbergeron/TreeshrewProject/iqtree/gtrees.tre"
genetree_names="/data/schwartzlab/zbergeron/TreeshrewProject/iqtree/gtrees.txt"
script_path="/data/schwartzlab/zbergeron/TreeshrewProject/screening/treescreen.R"

cd $SLURM_SUBMIT_DIR

module load R/4.0.3-foss-2020b

date
mkdir screening_assessments
Rscript ${script_path} ${concattree_path} ${genetree_path} ${genetree_names}
mv screening_output.csv screening_assessments/
date
