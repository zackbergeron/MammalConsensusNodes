#!/bin/bash
#SBATCH --job-name="amas"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=12   # processor core(s) per node
#SBATCH --mail-user="aknyshov@uri.edu"
#SBATCH --mail-type=END,FAIL

module purge
module load Python/3.7.4-GCCcore-8.3.0
path_to_amas="/home/aknyshov/alex_data/andromeda_tools/AMAS/amas/AMAS.py"
folder_with_loci="aligned_loci/"
cores_to_use=12

cd $SLURM_SUBMIT_DIR
date

mkdir amas_assessments

python ~/alex_data/tree_shew_analysis/TreeshrewProject/amas/run_amas.py ${folder_with_loci} ${cores_to_use} ${path_to_amas}

mv amas_total_results.txt amas_assessments/
rm amas_output_temp.txt

date