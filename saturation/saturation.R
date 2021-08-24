library(ape)
args = commandArgs(trailingOnly=TRUE)
argsLen <- length(args);
if (argsLen == 0){
  cat("Syntax: Rscript saturation.R [path to loci folder]\n")
  cat("Example: Rscript saturation.R ./dna/\n")
  quit()
}

dnaseqpath <- args[1]
print("get file list")
dnafiles <- list.files(path=dnaseqpath, pattern="\\.fasta$", full.names=FALSE, recursive=FALSE)
print("start analysis")
for (i in 1:length(dnafiles)){
	#path
	dnafile <- paste0(dnaseqpath,"/",dnafiles[i])
	#read
	alignment <- read.dna(dnafile, format="fasta")
	rawdist <- dist.dna(alignment, model = "raw", pairwise.deletion = T, as.matrix = T)
	tn93dist <- dist.dna(alignment, model = "TN93", pairwise.deletion = T, as.matrix = T)
	tn93dist2 <- tn93dist[lower.tri(tn93dist)]
	# get mean tn93 dist
	mean_dist <- mean(tn93dist2)
	sd_dist <- sd(tn93dist2)
	rawdist2 <- rawdist[lower.tri(rawdist)]
	regress <- lm(rawdist2 ~ tn93dist2)
	# get slope
	slope <- coef(regress)[2]
	# get r-squared
	Rsquared <- summary(regress)$r.squared
	output <- paste(dnafiles[i],slope,Rsquared,mean_dist,sd_dist,sep=",")
	write(x=output,file="saturation_output.csv", append=T)
}