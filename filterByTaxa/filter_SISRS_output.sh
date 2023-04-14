#!/bin/sh
#SBATCH --job-name="filter"
#SBATCH --time=30:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # processor core(s) per node - just use 1 - it do$
#SBATCH --mail-user="zbergeron@uri.edu" #CHANGE to your user email address
#SBATCH --mail-type=END,FAIL
### adjust/add sbatch flags as needed

cd $SLURM_SUBMIT_DIR

module purge
#for URI's Andromeda cluster
module load Biopython/1.78-foss-2020b

# adjust the command as needed, check the correct path to the python script
mkdir /data/schwartzlab/zbergeron/SISRS_mammals/filteredMammalLoci
python /data/schwartzlab/zbergeron/TreeshrewProject/filterByTaxa/filter_SISRS_output.py /data/schwartzlab/zbergeron/SISRS_mammals/Mammal_filterbytaxa_table.csv /data/schwartzlab/zbergeron/SISRS_mammals/SISRS_Run/aligned_contigs /data/schwartzlab/zbergeron/SISRS_mammals/filteredMammalLoci 0.33 18 4

date
