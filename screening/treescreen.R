library(ape)
library(phangorn)
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args);
if (argsLen == 0){
  cat("Syntax: Rscript treescreen.R [reference tree file] [path to loci folder]\n")
  cat("Example: Rscript treescreen.R tree.tre ./dna/\n")
  quit()
}

mltreepath <- args[1]
gtreespath <- args[2]
gtreenames <- args[3]

print("process concat trees")
concat_trees <- read.tree(mltreepath)
concat_trees_data <- lapply(concat_trees, function(tree) {
  return(setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which (y==x),y=tree$edge[,2])],tree$tip.label)  )
})
names(concat_trees_data) <- paste0("t",1:length(concat_trees))
concat_trees_data <- as.data.frame(concat_trees_data)
concat_trees_meanBL <- apply(concat_trees_data,1,mean)

write(x="locname,slope,Rsq",file="screening_output.csv")

gene_trees <- read.tree(gtreespath)
gene_names <- readLines(gtreenames)

print("process gene trees")
for (i in 1:length(gene_trees)){
	print (gene_names[i])
	gtree <- gene_trees[[i]]
	gtreeBL <- setNames(gtree$edge.length[sapply(1:length(gtree$tip.label),function(x,y) which (y==x),y=gtree$edge[,2])],gtree$tip.label)
	combolengths <- cbind(concat_trees_meanBL, gtreeBL = gtreeBL[names(concat_trees_meanBL)])
	regress <- lm(combolengths[,1] ~ combolengths[,2])
	# get slope
	slope <- coef(regress)[2]
	# get r-squared
	Rsquared <- summary(regress)$r.squared

	output <- paste(gene_names[i],slope,Rsquared,sep=",")
	write(x=output,file="screening_output.csv", append=T)
}