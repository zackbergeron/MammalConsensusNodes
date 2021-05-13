#!/bin/bash
#SBATCH --job-name="PMprep"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G

script_work_dir="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/"
array_script="/home/aknyshov/alex_data/tree_shew_analysis/TreeshrewProject/phylomad/phylomad_array.sh"
aligned_loci_path="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/aligned_loci/"
cd ${script_work_dir}
mkdir phylomad_assessment
cd phylomad_assessment
ls ${aligned_loci_path} | rev | cut -f1 -d/ | rev | split - aligned_loci_list_
arrayN=$(ls aligned_loci_list_* | wc -l)
ls aligned_loci_list_* > array_list.txt
echo sbatch --array=1-${arrayN}%40 ${array_script}