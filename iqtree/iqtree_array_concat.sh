#!/bin/bash
#SBATCH --job-name="IQcon"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH --mail-user="zbergeron@uri.edu"
#SBATCH --mail-type=ALL
#SBATCH -c 10
#SBATCH --mem-per-cpu=6G
#SBATCH --array=[1-102]%40

cd $SLURM_SUBMIT_DIR

date

module purge
module load Python/3.7.4-GCCcore-8.3.0

# first need to copy
# cp ../iqtree_phyloinference/array_list.txt .
# and
# cp ../iqtree_phyloinference/aligned_loci_list_* .

iqtree_assessment_path="/data/schwartzlab/zbergeron/SISRS_mammals/iqtree_assessment_mammals"
cp ${iqtree_assessment_path}/array_list.txt .
cp ${iqtree_assessment_path}/aligned_loci_list_* .

fileline=$(sed -n ${SLURM_ARRAY_TASK_ID}p array_list.txt)
aligned_loci_path="/data/schwartzlab/zbergeron/SISRS_mammals/filteredMammalLoci/"
iqtree_exe="/home/aknyshov/alex_data/andromeda_tools/iqtree-2.1.2-Linux/bin/iqtree2"
trees_to_eval="/data/schwartzlab/zbergeron/SISRS_mammals/hypothesis_trees/AllTrees.tree"
scripts_dir="/data/schwartzlab/zbergeron/TreeshrewProject/iqtree/"
path_to_amas="/data/schwartzlab/Biancani/AMAS/amas/AMAS.py"


infiles=$(cat ${fileline} | while read line; do echo ${aligned_loci_path}/${line}; done | paste -sd" ")

#amas concatenated
python3 ${path_to_amas} concat -f fasta -d dna --out-format fasta --part-format raxml -i $infiles -t concatenated_${SLURM_ARRAY_TASK_ID}.fasta -p partitions_${SLURM_ARRAY_TASK_ID}.txt


module purge
module load R/4.0.3-foss-2020b


Rscript ${scripts_dir}trimTrees.R concatenated_${SLURM_ARRAY_TASK_ID}.fasta ${trees_to_eval} ./trees_${SLURM_ARRAY_TASK_ID}.tre

sed -n 1p ./trees_${SLURM_ARRAY_TASK_ID}.tre > ./tree1_${SLURM_ARRAY_TASK_ID}.tre
${iqtree_exe} -nt 10 -s concatenated_${SLURM_ARRAY_TASK_ID}.fasta -spp partitions_${SLURM_ARRAY_TASK_ID}.txt -pre calcLnL1_${SLURM_ARRAY_TASK_ID} -g ./tree1_${SLURM_ARRAY_TASK_ID}.tre -m GTR+G -redo -wsl

sed -n 2p ./trees_${SLURM_ARRAY_TASK_ID}.tre > ./tree2_${SLURM_ARRAY_TASK_ID}.tre
${iqtree_exe} -nt 10 -s concatenated_${SLURM_ARRAY_TASK_ID}.fasta -spp partitions_${SLURM_ARRAY_TASK_ID}.txt -pre calcLnL2_${SLURM_ARRAY_TASK_ID} -g ./tree2_${SLURM_ARRAY_TASK_ID}.tre -m GTR+G -redo -wsl

sed -n 3p ./trees_${SLURM_ARRAY_TASK_ID}.tre > ./tree3_${SLURM_ARRAY_TASK_ID}.tre
${iqtree_exe} -nt 10 -s concatenated_${SLURM_ARRAY_TASK_ID}.fasta -spp partitions_${SLURM_ARRAY_TASK_ID}.txt -pre calcLnL3_${SLURM_ARRAY_TASK_ID} -g ./tree3_${SLURM_ARRAY_TASK_ID}.tre -m GTR+G -redo -wsl

${iqtree_exe} -nt 10 -s concatenated_${SLURM_ARRAY_TASK_ID}.fasta -spp partitions_${SLURM_ARRAY_TASK_ID}.txt -pre inference_${SLURM_ARRAY_TASK_ID} -m GTR+G -bb 1000 -alrt 1000 -redo -wsl -wsr

date
