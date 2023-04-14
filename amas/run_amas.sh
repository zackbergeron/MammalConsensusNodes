#!/bin/bash
#SBATCH --job-name="amas"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # processor core(s) per node
#SBATCH --mail-user="zbergeron@uri.edu"
#SBATCH --mail-type=END,FAIL

module purge
module load Python/3.7.4-GCCcore-8.3.0
path_to_run_amas_py="/data/schwartzlab/zbergeron/TreeshrewProject/amas/run_amas.py"
path_to_amas="/data/schwartzlab/Biancani/AMAS/amas/AMAS.py"
folder_with_loci="/data/schwartzlab/zbergeron/SISRS_mammals/filteredMammalLoci"
cores_to_use=36

cd $SLURM_SUBMIT_DIR
date

mkdir amas_assessments

python ${path_to_run_amas_py} ${folder_with_loci} ${cores_to_use} ${path_to_amas}

mv amas_total_results.txt amas_assessments/
rm amas_output_temp.txt

date
