#!/bin/sh
#SBATCH --job-name="filter"
#SBATCH --time=30:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - just use 1 
#SBATCH --mail-user="biancani@uri.edu" #CHANGE to your email
#SBATCH --mail-type=ALL

## UPDATE PATHS:
# path to SISRS loci (aligned contigs) folder:
LOCI=/data/schwartzlab/Biancani/PLACENTAL/SISRS_out/SISRS_Run/aligned_contigs
# path to folder containing filterByTaxa scripts:
SCRIPTS=/data/schwartzlab/Biancani/PlacentalPolytomy/annotation
# path to taxon group table (csv):
TXNGROUPS=/data/schwartzlab/Biancani/PlacentalPolytomy/filterByTaxa/groups.csv
# path to output folder for filtered loci (will be created by script if necessary):
OUTPUT=/data/schwartzlab/Biancani/PlacentalPolytomy/output/SISRS_out_filtered

## UPDATE parameters:
SEQCOMPLETE=0.33 # taxon sequence completeness, (e.g. 0.33 is 33% non N)
MINTAXA=18 # minimum number taxa to be present, e.g. 25
MINGROUPS=4 # minimum number of groups to be present, e.g. 4

cd $SLURM_SUBMIT_DIR

module purge
#for URI's Andromeda cluster
module load Biopython/1.78-foss-2020b 

mkdir -p $OUTPUT
python $SCRIPTS/filter_SISRS_output.py $TXNGROUPS $LOCI $OUTPUT $SEQCOMPLETE $MINTAXA $MINGROUPS