#!/bin/bash
#SBATCH --job-name="IQarr"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 10
#SBATCH --mem-per-cpu=6G
#SBATCH --array=[1-28]%28
#SBATCH --mail-user="biancani@uri.edu"
#SBATCH --mail-type=ALL

## update array line above based on output of iqtree prep script

cd $SLURM_SUBMIT_DIR

date
module load R/4.0.3-foss-2020b

# first need to copy 
# cp ../iqtree_phyloinference/array_list.txt .
# and
# cp ../iqtree_phyloinference/aligned_loci_list_* .

fileline=$(sed -n ${SLURM_ARRAY_TASK_ID}p array_list.txt)
aligned_loci_path="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/aligned_loci/"
iqtree_exe="/home/aknyshov/alex_data/andromeda_tools/iqtree-2.1.2-Linux/bin/iqtree2"
# trees_to_eval="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/treeshrew_trees_to_fit.tre"
trees_to_eval="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/primates_trees_to_fit.tre"
scripts_dir="/home/aknyshov/alex_data/tree_shew_analysis/TreeshrewProject/iqtree/"


infiles=$(cat ${fileline} | while read line; do echo ${aligned_loci_path}/${line}; done | paste -sd" ")

python3 ~/alex_data/andromeda_tools/AMAS/amas/AMAS.py concat -f fasta -d dna --out-format fasta --part-format raxml -i $infiles -t concatenated_${SLURM_ARRAY_TASK_ID}.fasta -p partitions_${SLURM_ARRAY_TASK_ID}.txt

Rscript ${scripts_dir}trimTrees.R concatenated_${SLURM_ARRAY_TASK_ID}.fasta ${trees_to_eval} ./trees_${SLURM_ARRAY_TASK_ID}.tre

${iqtree_exe} -nt 10 -s concatenated_${SLURM_ARRAY_TASK_ID}.fasta -spp partitions_${SLURM_ARRAY_TASK_ID}.txt -z ./trees_${SLURM_ARRAY_TASK_ID}.tre -pre calcLnL_${SLURM_ARRAY_TASK_ID} -n 0 -m GTR+G -wsl
${iqtree_exe} -nt 10 -s concatenated_${SLURM_ARRAY_TASK_ID}.fasta -spp partitions_${SLURM_ARRAY_TASK_ID}.txt -pre inference_${SLURM_ARRAY_TASK_ID} -m GTR+G -bb 1000 -alrt 1000 -wsr

date