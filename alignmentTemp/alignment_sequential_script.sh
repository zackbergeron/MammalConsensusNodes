#!/bin/sh

#SBATCH --job-name="align3"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH -c 4

### adjust/add sbatch flags as needed

cd $SLURM_SUBMIT_DIR

module load MAFFT/7.475-gompi-2020b-with-extensions
mkdir aligned_loci2

for infile in new_filtered_loci/*.fasta
do
	mafft --auto --thread 4 ${infile} > aligned_loci2/$(basename ${infile})
done
