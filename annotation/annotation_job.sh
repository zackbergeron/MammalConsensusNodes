#!/bin/bash
#SBATCH --job-name="annotate"
#SBATCH --time=100:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # processor core(s) per node
#SBATCH --mail-user="biancani@uri.edu" # CHANGE THIS to your user email address
#SBATCH --mail-type=ALL

## UPDATE:
# number of tasks/processor cores per node:
TASKS=36
# location of annotation scripts:
SCRIPTS=/data/schwartzlab/Biancani/PlacentalPolytomy/annotation
# path to FILTERED SISRS loci (aligned contigs) folder:
LOCI=/data/schwartzlab/Biancani/PlacentalPolytomy/output/SISRS_out_filtered
# path to output folder for annotation (will be created by script if necessary):
OUTPUT=/data/schwartzlab/Biancani/PlacentalPolytomy/output/annotation

# name of study taxon closest to reference:
taxonName="Pan_troglodytes"
# path to FASTA format reference genome:
GENOME="/data/schwartzlab/Biancani/PlacentalPolytomy/output/reference_genome/GCF_002880755.1_Clint_PTRv2_genomic.fna"
# path to GFF format reference annotation file:
ANNOTATION="/data/schwartzlab/Biancani/PlacentalPolytomy/output/reference_genome/GCF_002880755.1_Clint_PTRv2_genomic.gff"
# select output mode (either counts of different annotation types, or length proportion of each annotations type)
## (`c`) count the number of each feature per locus, for ex. 1 CDS, 2 introns, etc.
## (`l`) compute proportion of length of each feature type per locus, for ex. 0.2 CDS, 0.8 introns
outputMode="l"

mkdir $OUTPUT
cd $OUTPUT

#Andromeda (URI's cluster) specific
module purge
module load Biopython/1.78-foss-2020b
#
python ${SCRIPTS}/annotation_getTaxContigs.py ${taxonName} ${LOCI}

#Andromeda (URI's cluster) specific
module purge
module load BLAST+
#

makeblastdb -in ${GENOME} -dbtype nucl

blastn -query ${taxonName}.fasta -db ${GENOME} -outfmt 6 -num_threads ${TASKS} > blast_results.blast

python3 ${SCRIPTS}/annotation_blast_parser.py blast_results.blast > full_table.bed

sort -k1,1 -k2,2n full_table.bed > full_table_sorted.bed

#Andromeda (URI's cluster) specific
module purge
module load BEDTools/2.27.1-foss-2018b
#

bedtools intersect -a full_table_sorted.bed -b ${ANNOTATION} -wa -wb > full_table_annotated.bed

python3 ${SCRIPTS}/annotation_bed2table.py full_table_annotated.bed ${outputMode} > annotations.csv

date