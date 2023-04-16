#!/bin/bash
#SBATCH --job-name="IQarr"
#SBATCH --time=192:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH --mail-user="zbergeron@uri.edu"
#SBATCH --mail-type=ALL
#SBATCH -c 1
#SBATCH --mem-per-cpu=6G
#SBATCH --array=[1-102]%40

#specify working dirctory (the folder created by prep script)
array_work_dir="/data/schwartzlab/zbergeron/SISRS_mammals/iqtree_assessment_mammals/"

cd $array_work_dir

date
module load R/4.0.3-foss-2020b

fileline=$(sed -n ${SLURM_ARRAY_TASK_ID}p array_list.txt)
aligned_loci_path="/data/schwartzlab/zbergeron/SISRS_mammals/filteredMammalLoci/"
#location of iqtree executable:
iqtree_exe="/data/schwartzlab/alex/andromeda_tools/iqtree-2.1.2-Linux/bin/iqtree2"
#path to iqtree scripts
scripts_dir="/data/schwartzlab/zbergeron/TreeshrewProject/iqtree/"
# path to file containing 3 alternative hypotheses trees:
trees_to_eval="/data/schwartzlab/zbergeron/SISRS_mammals/hypothesis_trees/AllTrees.tree"

# Specify taxon list for each hypothesis tree:
# Tree 1 (Monophyletic Afrotheria) = Afrotheria List
focal_tips1="/data/schwartzlab/Biancani/PLACENTAL/hypothesis_trees/ZB_tips_Afrotheria.txt"
# List 2 (Monophyletic Boreoeutheria) = Boreoeutheria List
focal_tips2="/data/schwartzlab/Biancani/PLACENTAL/hypothesis_trees/ZB_tips_Boreoeutheria.txt"
# List 3 (Monophyletic Xenarthra) = Xenarthra List
focal_tips3="/data/schwartzlab/Biancani/PLACENTAL/hypothesis_trees/ZB_tips_Xenarthra.txt"
# Outgroup taxa list:
outgroup_tips="/data/schwartzlab/Biancani/PLACENTAL/hypothesis_trees/ZB_tips_Outgroup.txt"

> LnLs_${SLURM_ARRAY_TASK_ID}.csv
cat ${fileline} | while read line
do
  	echo $line
        Rscript ${scripts_dir}trimTrees.R ${aligned_loci_path}/${line} ${trees_to_eval} ./trees_${line}.tre

        sed -n 1p ./trees_${line}.tre > ./tree1_${line}.tre
        ${iqtree_exe} -nt 1 -t ./tree1_${line}.tre -s ${aligned_loci_path}/${line} --scf 500 --prefix concord1_${line}
        echo "ID,sCF,sCF_N,sDF1,sDF1_N,sDF2,sDF2_N,sN,debug" > scf/scf_${line}
        Rscript ${scripts_dir}getSCF.R concord1_${line}.cf.branch concord1_${line}.cf.stat ${focal_tips1} ${outgroup_tips} >> scf/scf_${line}
        ${iqtree_exe} -nt 1 -s ${aligned_loci_path}/${line} -pre calcLnL_${line} -g ./tree1_${line}.tre -m GTR+G -redo
        echo $line","$(grep "BEST SCORE FOUND" calcLnL_${line}.log | cut -f2 -d: ) >> LnLs_${SLURM_ARRAY_TASK_ID}.csv
        rm concord1_${line}.log concord1_${line}.cf.branch concord1_${line}.cf.stat concord1_${line}.cf.tree concord1_${line}.cf.tree.nex
        rm ./tree1_${line}.tre

        sed -n 2p ./trees_${line}.tre > ./tree2_${line}.tre
        ${iqtree_exe} -nt 1 -t ./tree2_${line}.tre -s ${aligned_loci_path}/${line} --scf 500 --prefix concord2_${line}
        Rscript ${scripts_dir}getSCF.R concord2_${line}.cf.branch concord2_${line}.cf.stat ${focal_tips2} ${outgroup_tips} >> scf/scf_${line}
        ${iqtree_exe} -nt 1 -s ${aligned_loci_path}/${line} -pre calcLnL_${line} -g ./tree2_${line}.tre -m GTR+G -redo
        echo $line","$(grep "BEST SCORE FOUND" calcLnL_${line}.log | cut -f2 -d: ) >> LnLs_${SLURM_ARRAY_TASK_ID}.csv
        rm concord2_${line}.log concord2_${line}.cf.branch concord2_${line}.cf.stat concord2_${line}.cf.tree concord2_${line}.cf.tree.nex
        rm ./tree2_${line}.tre

        sed -n 3p ./trees_${line}.tre > ./tree3_${line}.tre
        ${iqtree_exe} -nt 1 -t ./tree3_${line}.tre -s ${aligned_loci_path}/${line} --scf 500 --prefix concord3_${line}
        Rscript ${scripts_dir}getSCF.R concord3_${line}.cf.branch concord3_${line}.cf.stat ${focal_tips3} ${outgroup_tips} >> scf/scf_${line}
        ${iqtree_exe} -nt 1 -s ${aligned_loci_path}/${line} -pre calcLnL_${line} -g ./tree3_${line}.tre -m GTR+G -redo
        echo $line","$(grep "BEST SCORE FOUND" calcLnL_${line}.log | cut -f2 -d: ) >> LnLs_${SLURM_ARRAY_TASK_ID}.csv
        rm concord3_${line}.log concord3_${line}.cf.branch concord3_${line}.cf.stat concord3_${line}.cf.tree concord3_${line}.cf.tree.nex
        rm ./tree3_${line}.tre

        ${iqtree_exe} -nt 1 -s ${aligned_loci_path}/${line} -pre calcLnL_${line} -m GTR+G -redo
        echo $line","$(grep "BEST SCORE FOUND" calcLnL_${line}.log | cut -f2 -d: ) >> LnLs_${SLURM_ARRAY_TASK_ID}.csv

        rm ./trees_${line}.tre calcLnL_${line}.ckp.gz calcLnL_${line}.iqtree calcLnL_${line}.log calcLnL_${line}.treefile calcLnL_${line}.trees
        rm calcLnL_${line}.uniqueseq.phy
done
