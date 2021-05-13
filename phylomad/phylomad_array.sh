#!/bin/bash
#SBATCH --job-name="PMarr"
#SBATCH --time=72:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G


module load R/4.0.3-foss-2020b
date

fileline=$(sed -n ${SLURM_ARRAY_TASK_ID}p array_list.txt)
aligned_loci_path="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing/aligned_loci/"
base_phylomad_script="/home/aknyshov/alex_data/andromeda_tools/phylomad/codeFolder/commandLineFile.Rscript"

cat ${fileline} | while read line
do
	echo $line #locus file
	sed -e 's+input$dataPath <- c(".*")+input$dataPath <- c("'"$aligned_loci_path""$line"'")+g' ${base_phylomad_script} | sed -e 's/alNames <- c(".*")/alNames <- c("'"$line"'")/g' | sed -e 's+input$outputFolder <- ".*"+input$outputFolder <- "'"$PWD"'/"+g' > script_${SLURM_ARRAY_TASK_ID}.R
	Rscript script_${SLURM_ARRAY_TASK_ID}.R
	rm script_${SLURM_ARRAY_TASK_ID}.R
done