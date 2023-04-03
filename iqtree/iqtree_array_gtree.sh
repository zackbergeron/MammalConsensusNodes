#!/bin/bash
#SBATCH --job-name="IQarr"
#SBATCH --time=48:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G
#SBATCH --array=[1-28]%28
#SBATCH --mail-user="biancani@uri.edu"
#SBATCH --mail-type=ALL

## update array line above based on output of iqtree prep script


cd $SLURM_SUBMIT_DIR

date

fileline=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../iqtree_assessment/array_list.txt)
aligned_loci_path="/home/aknyshov/alex_data/tree_shew_analysis/SISRS/post_processing_v2/aligned_loci/"
iqtree_exe="/home/aknyshov/alex_data/andromeda_tools/iqtree-2.1.2-Linux/bin/iqtree2"

cat ../iqtree_assessment/${fileline} | while read line
do
	echo $line
	${iqtree_exe} -nt 1 -s ${aligned_loci_path}/${line} -pre inference_${line} -alrt 1000 -m GTR+G
	rm inference_${line}.ckp.gz inference_${line}.iqtree inference_${line}.log inference_${line}.bionj inference_${line}.mldist inference_${line}.uniqueseq.phy
done
