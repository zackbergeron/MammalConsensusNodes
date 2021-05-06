#!/bin/sh

#SBATCH --job-name="filter"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)

### adjust/add sbatch flags as needed

cd $SLURM_SUBMIT_DIR

module purge
#for URI's Andromeda cluster
module load Biopython/1.78-foss-2020b 

# adjust the command as needed, check the correct path to the python script
mkdir <path for output folder>
python <path to>filter_SISRS_output.py <path to taxa group table> <path to SISRS loci folder> <path for output folder> <taxon sequence completeness, e.g. 0.33> <num taxa to be present, e.g. 25> <number of groups to be present, e.g. 6>