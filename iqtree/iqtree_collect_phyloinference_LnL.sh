#!/bin/bash
#SBATCH --job-name="iqout"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH --mail-user="zbergeron@uri.edu"
#SBATCH --mail-type=ALL
#SBATCH -c 1
#SBATCH --mem-per-cpu=8G
#SBATCH --array=[1-102]%40

cd $SLURM_SUBMIT_DIR

arrayLen=$(ls partitions_*.txt|wc -l)
assessment_folder="/data/schwartzlab/zbergeron/TreeshrewProject/iqtree"

> combined_iqtree_dLnLs_concat.csv
date
for f in $(seq 1 ${arrayLen})
do
	cat ${assessment_folder}/partitions_${f}.txt | while read l
	do
		locname=$(echo ${l} | cut -f2 -d" " | cut -f2- -d_)
		range1=$(echo ${l} | cut -f4 -d" ")
		tree1=$(sed -n 2p ${assessment_folder}/calcLnL1_${f}.sitelh | awk -v a="${range1}" 'BEGIN {split(a, A, /-/)} {x=0;for(i=A[1]+1;i<=A[2]+1;i++)x=x+$i;print x}')
		tree2=$(sed -n 2p ${assessment_folder}/calcLnL2_${f}.sitelh | awk -v a="${range1}" 'BEGIN {split(a, A, /-/)} {x=0;for(i=A[1]+1;i<=A[2]+1;i++)x=x+$i;print x}')
		tree3=$(sed -n 2p ${assessment_folder}/calcLnL3_${f}.sitelh | awk -v a="${range1}" 'BEGIN {split(a, A, /-/)} {x=0;for(i=A[1]+1;i<=A[2]+1;i++)x=x+$i;print x}')
		tree4=$(sed -n 2p ${assessment_folder}/inference_${f}.sitelh | awk -v a="${range1}" 'BEGIN {split(a, A, /-/)} {x=0;for(i=A[1]+1;i<=A[2]+1;i++)x=x+$i;print x}')
		echo ${locname},${tree1},${tree2},${tree3},${tree4} >> combined_iqtree_dLnLs_concat.csv
	done
done
date
