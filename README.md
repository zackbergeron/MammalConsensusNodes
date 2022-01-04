# TreeshrewProject
Scripts associated with the Treeshrew project

## Filter by taxa
This step can be done after the next one, for example when loci are already aligned as part of the SISRS pipeline
See **filterByTaxa** folder
* filter_SISRS_output.sh - slurm script to run
* filter_SISRS_output.py - script run by the previous shell script
Loci are first filtered from incomplete sequences (for the treeshrew project each passing sequence must have over 33% non N), then are filtered by the number of taxa retained (for the treeshrew project the locus passes if 25 taxa are remaining) and then by the number of groups present (for the treeshrew project all 6 groups must be present).
The taxa to group correspondence table is expected to be a csv file and look like this:
```
group1,taxonName1
group1,taxonName2
group2,taxonName3
group2,taxonName4
```


## Align loci
Obsolete scripts - SISRS now includes a step to align obtained loci
See **alignmentTemp** folder
* alignment_prep_script.sh - to prep the slurm array job
* alignment_array_script.sh - to run the alignment array job

## Locus annotation

This step aligns SISRS loci of a particular taxon to the reference sequence (ideally, of the same taxon). A custom script is used to retreive a reference taxon for the loci. Then BLAST is run. A custom python script is used to filter the output and convert it to BED. Overlapping hits of similar scores as well as very disjunct alignments are discarded. The BED file is then sorted and intersected with the GFF file for the reference sequence. A custom python script then processess the intersected BED file to produce the final output, either counts of different types of annotations, or length proportions of each of types of annotations. Currently the following types of annotations are recorded: pseudogene, CDS, UTR, intron, lnc_RNA, other (any other type), unannotated (or intergenic). See **annotation** folder
* annotation_job.sh - slurm script to submit; adjust file paths; adjust the command options on line 30 for the last script (see below)
* annotation_blast_parser.py - BLAST to BED script run by the previous shell script; adjust parameters on lines 10-12 as needed.
* annotation_bed2table.py - BED to table script run by the previous shell script; two options are provided: count the number of each feature per locus (`c`), for ex. 1 CDS, 2 introns, etc.; or compute proportion of length of each feature type per locus (`l`), for ex. 0.2 CDS, 0.8 introns.
* getTaxContigs.py - script to extract SISRS loci of a particular taxon into a single file for use with BLAST

## Assess locus properties

### AMAS
This step runs AMAS to assess locus features. Since AMAS takes file names as command line arguments, there's a limit on how long the line can be, and we have several hundred thousand files, AMAS is run in batches by the driver script. See **amas** folder
* run_amas.sh - slurm script to run; adjust `12` to the number of cores you wish to use
* run_amas.py - script run by the previous shell script

### Saturation
This step assesses uncorrected (raw) and corrected (TN93) pairwise distances for each locus, fits linear model to the data, and outputs the slope and R squared (modified from Borowiec et al. 2015). If the fit is good, the slope would be indicative of the degree of saturation. See **saturation** folder
* saturation.sh - slurm script to run; adjust paths to loci and to the saturation.R script.
* saturation.R - script run by the previous shell script

### Phylomad
This step runs phylomad. As phylomad has "control file" that specifies all parameters of the run, it is problematic to adjust it for each of the several hundred thousand runs needed. The solution here is to first adjust all parameters in the control file besides the filenames. Then a driver script takes care of adjusting the control file for each locus and running phylomad on it. Because phylomad is computationally intensive and we have too many alignments to run through, this step is set up as an array job. See **phylomad** folder
* phylomad_prep.sh - slurm script to set up the phylomad array; adjust the paths, possibly also the number of simultaneously run jobs on line 18 after the % sign.
* phylomad_array.sh - slurm script to run the array job; adjust the paths.
* phylomad_collect_output.sh - slurm script to collect the results, expects results to be in the folder phylomad_assessment, adjust the script otherwise
* phylomad_collect_array.sh - in case of large dataset, use this array script to collect the results instead, and then concatenate array results.

### Taxon composition
This step counts number of taxa in each group as defined for the filter by taxa step. See **misc**
* taxon_composition_slurm.sh - slurm script to run

## Assess the phylogenetic signal

### Obtain the trees to score

### Run IQ-TREE
This step runs several assessments using IQ-TREE and several custom scripts. See **iqtree** folder
* iqtree_prep.sh - slurm script to run: sets up the folders, lists of files to process and produces the command to run after; adjust paths on lines 9-12 as well as the number of simultaneous jobs on line 20 as needed.
* iqtree_array.sh - slurm script to submit: runs the analyses; **submit using the command output by the previous step in its log file (\*.out)**; adjust paths on lines 16-27 as needed.
* trimTrees.R - script to trim taxa off of the main tree(s) depending on the taxon composition of a particular alignment, run by the previous shell script.
* getSCF.R - script to extract sCF values for the branch of interest from IQ-TREE output, run by the previous shell script.

