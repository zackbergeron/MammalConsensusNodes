#!/bin/bash
#SBATCH --job-name="annotate"
#SBATCH --time=48:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node
#SBATCH --mail-user="aknyshov@uri.edu"
#SBATCH --mail-type=END,FAIL

module purge
module load BLAST+

cd $SLURM_SUBMIT_DIR
date

fileToMap="../refs_for_annotation/Homo_sapiens.fas"
assemblyDB="/home/aknyshov/alex_data/tree_shew_analysis/NCBI_genomes/Homo_sapiens/GCF_000001405.39_GRCh38.p13_genomic.mod.fna"
assemblyGFF="/home/aknyshov/alex_data/tree_shew_analysis/NCBI_genomes/Homo_sapiens/GCF_000001405.39_GRCh38.p13_genomic.gff"

blastn -query ${fileToMap} -db ${assemblyDB} -outfmt 6 -num_threads 20 > blast_results.blast

python3 ~/alex_data/tree_shew_analysis/TreeshrewProject/annotation/annotation_blast_parser.py blast_results.blast > full_table.txt

sort -k1,1 -k2,2n full_table.txt > full_table_sorted.bed

module purge
module load BEDTools/2.27.1-foss-2018b

bedtools intersect -a full_table_sorted.bed -b ${assemblyGFF} -wa -wb > full_table_annotated.bed

python3 ~/alex_data/tree_shew_analysis/TreeshrewProject/annotation/annotation_bed2table2.py full_table_annotated.bed > annotations.txt

date
