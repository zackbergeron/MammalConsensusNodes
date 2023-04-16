#!/bin/bash
#SBATCH --job-name="IQout"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH --mail-user="zbergeron@uri.edu"
#SBATCH --mail-type=ALL
#SBATCH -c 1
#SBATCH --mem-per-cpu=8G

cd $SLURM_SUBMIT_DIR

# submit after iqtree_array_gtree.sh script

date
> gtrees.txt; cat ../../SISRS_mammals/iqtree_assessment_mammals/array_list.txt | while read line1; do cat ../../SISRS_mammals/iqtree_assessment_mammals/${line1} >> gtrees.txt; done
> gtrees.tre; cat gtrees.txt | while read line; do cat ./inference_${line}.treefile >> gtrees.tre; done
date
