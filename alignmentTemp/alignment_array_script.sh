#!/bin/sh

#SBATCH --job-name="align2"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)
#SBATCH -c 4

### adjust/add sbatch flags as needed

cd $SLURM_SUBMIT_DIR

module load MAFFT/7.475-gompi-2020b-with-extensions

fileline=$(sed -n ${SLURM_ARRAY_TASK_ID}p align_list.txt)
cat alignmentGroups/${fileline} | while read line
do
	mafft --auto --thread 4 new_filtered_loci/$line > aligned_loci/$line
done
