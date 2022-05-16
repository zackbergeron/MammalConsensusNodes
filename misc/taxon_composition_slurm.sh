#!/bin/bash
#SBATCH --job-name="taxastat"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G

cd $SLURM_SUBMIT_DIR
alignments_folder="aligned_loci"
taxon_table="../../cactus_mammal_alignment/taxalist.csv"
echo "Alignment_name,"$(cut -f1 -d, ${taxon_table} | uniq -c | awk '{print $2}' | paste -sd",") > taxon_composition.csv
for f in ${alignments_folder}/*.fasta
do
	echo $(echo ${f} | cut -f2 -d/)","$(grep ">" ${f} | sed -e 's/>//g' | grep -f - ${taxon_table} | cut -f1 -d, | uniq -c | awk '{print $1}' | paste -sd",") >> taxon_composition.csv
done