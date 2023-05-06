# Mammal Consensus Nodes
Scripts associated with the [Mammal Consensus Node Project](https://github.com/zackbergeron/MammalConsensusNodes).
Adapted from [TreeshrewProject](https://github.com/AlexKnyshov/TreeshrewProject) and [PlacentalPolytomyProject](https://github.com/LMBiancani/PlacentalPolytomy).

Scripts associated with the Placental Mammal Project

Filter by taxa

This step can be done after loci are already aligned as part of the SISRS pipeline See filterByTaxa folder

filter_SISRS_output.sh - slurm script to run
filter_SISRS_output.py - script run by the previous shell script Loci are first filtered from incomplete sequences (for the treeshrew project each passing sequence must have over 33% non N), then are filtered by the number of taxa retained (for the treeshrew project the locus passes if 25 taxa are remaining) and then by the number of groups present (for the treeshrew project all 6 groups must be present). The taxa to group correspondence table is expected to be a csv file and look like this:
group1,taxonName1
group1,taxonName2
group2,taxonName3
group2,taxonName4
