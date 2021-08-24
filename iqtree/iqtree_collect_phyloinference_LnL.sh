#!/bin/bash
#SBATCH --job-name="iqout"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=8G

cd $SLURM_SUBMIT_DIR

arrayLen=387
assessment_folder="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/iqtree_phyloinference_treeshrew"

> combined_iqtree_dLnLs_concat.csv
date
for f in $(seq 1 ${arrayLen})
do
	cat ${assessment_folder}/partitions_${f}.txt | while read l
	do
		locname=$(echo ${l} | cut -f2 -d" " | cut -f2- -d_)
		range1=$(echo ${l} | cut -f4 -d" ")
		tree1=$(sed -n 2p ${assessment_folder}/calcLnL_${f}.sitelh | awk -v a="${range1}" 'BEGIN {split(a, A, /-/)} {x=0;for(i=A[1]+1;i<=A[2]+1;i++)x=x+$i;print x}')
		tree2=$(sed -n 3p ${assessment_folder}/calcLnL_${f}.sitelh | awk -v a="${range1}" 'BEGIN {split(a, A, /-/)} {x=0;for(i=A[1]+1;i<=A[2]+1;i++)x=x+$i;print x}')
		tree3=$(sed -n 4p ${assessment_folder}/calcLnL_${f}.sitelh | awk -v a="${range1}" 'BEGIN {split(a, A, /-/)} {x=0;for(i=A[1]+1;i<=A[2]+1;i++)x=x+$i;print x}')
		echo ${locname},${tree1},${tree2},${tree3} >> combined_iqtree_dLnLs_concat.csv
	done
done
date