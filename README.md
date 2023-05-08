# Mammal Consensus Nodes
Scripts associated with the [Mammal Consensus Node Project](https://github.com/zackbergeron/MammalConsensusNodes).
Adapted from [TreeshrewProject](https://github.com/AlexKnyshov/TreeshrewProject) and associated with the [PlacentalPolytomyProject](https://github.com/LMBiancani/PlacentalPolytomy).

Scripts associated with the Placental Mammal Project

Filter by taxa

This step can be done after loci are already aligned as part of the SISRS pipeline. See filterByTaxa folder

filter_SISRS_output.sh - slurm script to run
filter_SISRS_output.py - script run by the previous shell script Loci are first filtered from incomplete sequences (for the treeshrew project each passing sequence must have over 33% non N), then are filtered by the number of taxa retained (for the treeshrew project the locus passes if 25 taxa are remaining) and then by the number of groups present (for the treeshrew project all 6 groups must be present). The taxa to group correspondence table is expected to be a csv file and look like this:
group1,taxonName1
group1,taxonName2
group2,taxonName3
group2,taxonName4



Locus filtering using Branch Length Correlation

This step can be used to filter out outlier loci based on discordance in branch length distribution. For use with this script a gene tree for each locus needs to be computed. See IQ-TREE analyses below, and iqtree_array_gtree.sh and iqtree_collect_gtrees.sh in particular. Output is a table with regression slope and R-squared, the latter can be used to rank and filter out loci with worst (lowest) values. See screening folder

treescreen.sh - to run the screening job
treescreen.R - tree screening R script, run by the shell script above



Assess locus properties

AMAS

This step runs AMAS to assess locus features. Since AMAS takes file names as command line arguments, there's a limit on how long the line can be, and we have several hundred thousand files, AMAS is run in batches by the driver script. See amas folder

run_amas.sh - slurm script to run; adjust 12 to the number of cores you wish to use
run_amas.py - script run by the previous shell script



Assess the phylogenetic signal

Obtain the trees to score

Run IQ-TREE

This step runs several assessments using IQ-TREE and several custom scripts. See iqtree folder

iqtree_prep.sh - slurm script to run: sets up the folders, lists of files to process and produces the command to run after; adjust paths on lines 9-12 as well as the number of simultaneous jobs on line 20 as needed.
iqtree_array.sh - slurm script to submit: runs the analyses; submit using the command output by the previous step in its log file (*.out); adjust paths on lines 16-27 as needed.
iqtree_array_concat.sh - slurm script to submit: infer concatenation trees
iqtree_array_gtree.sh - slurm script to submit: infer gene trees
iqtree_collect_gtrees.sh - slurm script to submit: collect gene tree data
iqtree_collect_output.sh - slurm script to submit: collect individual fit assessment data
iqtree_collect_phyloinference_LnL.sh - slurm script to submit: collect concatenation fit assessment data
trimTrees.R - script to trim taxa off of the main tree(s) depending on the taxon composition of a particular alignment, run by the previous shell script.
getSCF.R - script to extract sCF values for the branch of interest from IQ-TREE output, run by the previous shell script.
