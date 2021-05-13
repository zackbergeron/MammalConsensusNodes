library (ape)
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args);
if (argsLen == 0){
  cat("Syntax: Rscript trimTrees.R [path to dna alignment] [path to trees file] [output tree path]\n")
  cat("Example: Rscript trimTrees.R dna.fas trees.tre trees_trimmed.tre\n")
  quit()
}
alignment <- read.dna(args[1], format="fasta")
trees <- read.tree(args[2])
outtreesname <- args[3]
taxapresent <- rownames(alignment)
outtrees <- list()
i <- 1
for (t in trees){
	treetips <- t$tip.label
	tipstoremove <- treetips[!(treetips %in% taxapresent)]
	outtrees[[i]] <- drop.tip(t, tipstoremove)
    i <- i + 1
}
class(outtrees) <- "multiPhylo"
write.tree(outtrees, outtreesname)