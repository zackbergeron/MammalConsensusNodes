#!/bin/bash
#SBATCH --job-name="amas"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # processor core(s) per node
#SBATCH --mail-user="biancani@uri.edu"
#SBATCH --mail-type=ALL

# UPDATE:
# path to AMAS program (AMAS.py):
AMAS=/data/schwartzlab/Biancani/AMAS/amas/AMAS.py
# number of tasks/cores per node:
CORES=36
# location of amas scripts:
SCRIPTS=/data/schwartzlab/Biancani/PlacentalPolytomy/amas
# path to FILTERED SISRS loci (aligned contigs) folder:
LOCI=/data/schwartzlab/Biancani/PlacentalPolytomy/output/SISRS_out_filtered
# path to output folder for amas (will be created by script if necessary):
OUTPUT=/data/schwartzlab/Biancani/PlacentalPolytomy/output/amas

#Andromeda (URI's cluster) specific
module purge
module load Python/3.7.4-GCCcore-8.3.0
#

mkdir -p ${OUTPUT}/amas_assessments
cd ${OUTPUT}
date

python ${SCRIPTS}/run_amas.py ${LOCI} ${CORES} ${AMAS}

mv amas_total_results.txt amas_assessments/
rm amas_output_temp.txt

date