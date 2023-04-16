#!/bin/bash
#SBATCH --job-name="IQarr"
#SBATCH --time=48:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH --mail-user="zbergeron@uri.edu"
#SBATCH --mail-type=ALL
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G
#SBATCH  --array=[1-102]%40

module load IQ-TREE/2.1.2-foss-2020a
cd $SLURM_SUBMIT_DIR

# make sure to set up as array job with last flag above

date

fileline=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../../SISRS_mammals/iqtree_assessment_mammals/array_list.txt)
aligned_loci_path="/data/schwartzlab/zbergeron/SISRS_mammals/filteredMammalLoci/"

cat ../../SISRS_mammals/iqtree_assessment_mammals/${fileline} | while read line
do
  	echo $line
        iqtree2 -nt 1 -s ${aligned_loci_path}/${line} -pre inference_${line} -alrt 1000 -m GTR+G
        rm inference_${line}.ckp.gz inference_${line}.iqtree inference_${line}.log inference_${line}.bionj inference_${line}.mldist inference_${line}.uniqueseq.phy
done
