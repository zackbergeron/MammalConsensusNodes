library (ape)
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args);
if (argsLen == 0){
  cat("Syntax: Rscript getSCF.R [path to .cf.branch] [path to .cf.stat] [path to tips list] [path to outgroup tips list] \n")
  cat("Example: Rscript getSCF.R test.cf.branch test.cf.stat tipsToAnalyse.txt outgroup.txt\n")
  quit()
}
reftree <- read.tree(args[1])
stattable <- read.table(args[2],header=TRUE, fill=TRUE)
filtertaxa <- readLines(args[3])
outtaxa <- readLines(args[4])

#since IQ-TREE occasionally reroots trees, have to reroot back by outgroup in order for the results to make sense

presentouttaxa <- reftree$tip.label[reftree$tip.label %in% outtaxa]
if (length(presentouttaxa) == 1){
	outgroupTipID <- which(reftree$tip.label == presentouttaxa)
	outgroupMRCA <- reftree$edge[reftree$edge[,2] == outgroupTipID,1]
} else {
	outgroupMRCA <- getMRCA(reftree, presentouttaxa)
}
if (outgroupMRCA != reftree$edge[1,1]){
reftree <- root(reftree, node= outgroupMRCA, edgelabel = T, resolve.root=TRUE) #---> this was change needed to root with our trees
}

presentTestTaxa <- reftree$tip.label[reftree$tip.label %in% filtertaxa]
if (length(presentTestTaxa) <= 1) {
	SCFofInterest <- c(" ", " ", " ", " ", " ", " ", " ", " ")
	debug_opt <- 1
} else {
	taxonNode <- getMRCA(reftree, presentTestTaxa)
	nodeToRetrieve<- as.numeric(reftree$node.label[taxonNode-length(reftree$tip.label)])
	debug_opt <- 2
	SCFofInterest <- stattable[stattable$ID == nodeToRetrieve,1:8]
}

write(paste(c(SCFofInterest,debug_opt),collapse=","),file="")
