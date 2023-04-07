#!/bin/bash
#SBATCH --job-name="IQprep"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-user="biancani@uri.edu"
#SBATCH --mail-type=ALL

# UPDATE:
# path to main output directory:
script_work_dir=/data/schwartzlab/Biancani/PlacentalPolytomy/output
# name output folder for iqtree (will be created by script if necessary):
array_work_folder="iqtree_assessment"
# location of iqtree scripts:
scripts_dir=/data/schwartzlab/Biancani/PlacentalPolytomy/iqtree
# path to FILTERED SISRS loci (aligned contigs) folder:
aligned_loci_path=/data/schwartzlab/Biancani/PlacentalPolytomy/output/SISRS_out_filtered
# number of simultaneous tasks for subsequent array jobs:
array_tasks=40

cd ${script_work_dir}
mkdir -p ${array_work_folder}
cd ${array_work_folder}
mkdir scf
# extract filenames from aligned_loci_path and split into bins of 4000
ls ${aligned_loci_path} | rev | cut -f1 -d/ | rev | split -l 4000 - aligned_loci_list_
arrayN=$(ls aligned_loci_list_* | wc -l)
ls aligned_loci_list_* > array_list.txt
if [ $arrayN -lt $array_tasks ]
    then
      array_tasks=$arrayN
fi

echo "#SBATCH --array=[1-${arrayN}]%${array_tasks}"
