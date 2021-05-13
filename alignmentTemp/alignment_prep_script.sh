#!/bin/sh

#SBATCH --job-name="align1"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)

### adjust/add sbatch flags as needed

cd $SLURM_SUBMIT_DIR

mkdir alignmentGroups
ls new_filtered_loci | split - align_group
mv align_group* alignmentGroups
ls alignmentGroups > align_list.txt
mkdir aligned_loci