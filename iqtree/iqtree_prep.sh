#!/bin/bash
#SBATCH --job-name="IQprep"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G

script_work_dir="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/"
scripts_dir="/home/aknyshov/alex_data/tree_shew_analysis/TreeshrewProject/iqtree/"
aligned_loci_path="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/aligned_loci/"
array_work_folder="iqtree_assessment_primates"
cd ${script_work_dir}
mkdir ${array_work_folder}
cd ${array_work_folder}
mkdir scf
ls ${aligned_loci_path} | rev | cut -f1 -d/ | rev | split - aligned_loci_list_
arrayN=$(ls aligned_loci_list_* | wc -l)
ls aligned_loci_list_* > array_list.txt
echo sbatch --array=1-${arrayN}%40 ${scripts_dir}iqtree_array.sh