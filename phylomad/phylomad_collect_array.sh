#!/bin/bash
#SBATCH --job-name="PMout"
#SBATCH --time=6:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G


date

array_list_line=$(sed -n ${SLURM_ARRAY_TASK_ID}p phylomad_assessment/array_list.txt)

> combined_phylomad_${array_list_line}.csv
cat phylomad_assessment/${array_list_line} | while read aligned_loci_list_line
do
	f="phylomad_assessment/"${aligned_loci_list_line}".phylomad.subst"
	loc=$(head -1 ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	val1=$(grep "chisq.stdev.from.pred.dist" ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	val2=$(grep "multlik.stdev.from.pred.dist" ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	val3=$(grep "biochemdiv.stdev.from.pred.dist" ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	val4=$(grep "consind.stdev.from.pred.dist" ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	val5=$(grep "delta.stdev.from.pred.dist" ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	val6=$(grep -w "brsup.stdev.from.pred.dist" ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	val7=$(grep "CIbrsup.stdev.from.pred.dist" ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	val8=$(grep "trlen.stdev.from.pred.dist" ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	val9=$(grep "maha.stdev.from.pred.dist" ${f}/output.pvals.PhyloMAd.csv | cut -f2 -d,)
	echo ${loc},${val1},${val2},${val3},${val4},${val5},${val6},${val7},${val8},${val9} >> combined_phylomad_${array_list_line}.csv
done