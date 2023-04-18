title: "Placental Mammal Project Data Analysis"
author: "Zachary Bergeron, credit to Alex Knyshov for majority of base script"
date: "03/31/2023"
output: html_document

  
library(tidyr)
library(dplyr)
library(scaleboot)

###Beginning Notes
  #'PRI' = Primates from Alex; here it corresponds to the 3 Monophyletic Placental Clades Constraint Trees + Star
  #'TSH' = Treeshreew from Alex; here it corresponds to the 3 alternative hypotheses for root of placental mammals
  #'i' after 'PRI' = individual fit
  #'c' after 'PRI' = concatenated fit
  #'r' after 'i/c' = rescaled (i.e. dividied by locus length)
  #'A' = Monophyletic Afrotheria constraint tree/Comparison w/ Star
  #'B' = Monophyletic Boreoeutheria constraint tree/Comparison w/ Star
  #'X' = Monophyletic Xenarthra constraint tree/Comparison w/ Star
  
###Individual fit data prep; Monophyletic clade constraints

dLnLtabPRIi <- read.csv("combined_iqtree_mammals_dLnLs.csv", header = F, stringsAsFactors = F)
#deleted last 9 LnL values as they were repeats from previous check runs; 405,629 to 405,620

dLnLtabPRIi <- dLnLtabPRIi[1:405620,]
dLnLtabPRIi$Tree <- rep(1:4,101405)

#had to go wide-format b/c original csv had long-format
dLnLtabPRIi_wide <- pivot_wider(dLnLtabPRIi, names_from = Tree, values_from = V2)

#didn't screen like Alex; skipped

#AU test to assess signficance of LnL difference across whole dataset for any one comparison 

#Mono-Afro vs. Star
library(scaleboot)
AUtestPRIiA <- relltest(na.omit(as.matrix(dLnLtabPRIi_wide[,c(2,5)])),nb=1000)
AUtestPRIisummaryA <- summary(AUtestPRIiA)
AUtestPRIisummaryA

AUtestPRIisummaryAttrA <- attributes(AUtestPRIisummaryA)
AUtestPRIipvalA <- 1-AUtestPRIisummaryAttrA$table$value[which.max(AUtestPRIisummaryAttrA$table$value[,"sk.3"]),"sk.3"]
AUtestPRIipvalA

#Mono-Boreo vs. Star
AUtestPRIiB <- relltest(na.omit(as.matrix(dLnLtabPRIi_wide[,c(3,5)])),nb=1000)
AUtestPRIisummaryB <- summary(AUtestPRIiB)
AUtestPRIisummaryB

AUtestPRIisummaryAttrB <- attributes(AUtestPRIisummaryB)
AUtestPRIipvalB <- 1-AUtestPRIisummaryAttrB$table$value[which.max(AUtestPRIisummaryAttrB$table$value[,"sk.3"]),"sk.3"]
AUtestPRIipvalB

#Mono-Xena vs. Star
AUtestPRIiX <- relltest(na.omit(as.matrix(dLnLtabPRIi_wide[,c(4,5)])),nb=1000)
AUtestPRIisummaryX <- summary(AUtestPRIiX)
AUtestPRIisummaryX

AUtestPRIisummaryAttrX <- attributes(AUtestPRIisummaryX)
AUtestPRIipvalX <- 1-AUtestPRIisummaryAttrX$table$value[which.max(AUtestPRIisummaryAttrX$table$value[,"sk.3"]),"sk.3"]
AUtestPRIipvalX

#rescaling (just dividing by locus length); added [,5] line b/c I have extra tree
amatab <- read.table("amas_total_results.txt", header = T, stringsAsFactors = F, check.names = F)
amatabfiltPRIi <- amatab[match (dLnLtabPRIi$V1,amatab$Alignment_name ),]
dLnLtabPRIir <- dLnLtabPRIi_wide
dLnLtabPRIir[,2] <- as.numeric(dLnLtabPRIir[,2])
dLnLtabPRIir[,3] <- as.numeric(dLnLtabPRIir[,3])
dLnLtabPRIir[,4] <- as.numeric(dLnLtabPRIir[,4])
dLnLtabPRIir[,5] <- as.numeric(dLnLtabPRIir[,5])
dLnLtabPRIir <- na.omit(dLnLtabPRIir)
dLnLtabPRIir[,2:5] <- (dLnLtabPRIir[,2:5] - apply(dLnLtabPRIir[,2:5], 1, min))/amatabfiltPRIi$Alignment_length
colnames(dLnLtabPRIir) <- c("locname","MonoA","MonoB","MonoX","Star")

#get outlier/nonoutlier data

#this was Alex's code but decided to simplify
#dLnLtabPRIirmax <- apply(dLnLtabPRIir[,2:5], 1, max)
#dLnLtabPRIirOUT <- dLnLtabPRIir[order(dLnLtabPRIirmax, decreasing = T),][1:1000,]

#can compare likelihoods b/c they were rescaled (i.e. divided by loci length)
#'outliers' considered "good" loci here as they provide strong signal for one of the hypotheses

#Afrotheria
dLnLtabPRIirOUTA <- dLnLtabPRIir %>% 
  rowwise() %>% mutate(max25 = max(MonoA,Star)) %>%
  arrange(desc(max25)) %>% head(1000) %>% mutate(dif = MonoA-Star)

dLnLtabPRIirNONOUTA <- dLnLtabPRIir[!(dLnLtabPRIir$locname %in% dLnLtabPRIirOUTA$locname),]

#Boreoeutheria
dLnLtabPRIirOUTB <- dLnLtabPRIir %>% 
  rowwise() %>% mutate(max35 = max(MonoB,Star)) %>%
  arrange(desc(max35)) %>% head(1000) %>% mutate(dif = MonoB-Star)

dLnLtabPRIirNONOUTB <- dLnLtabPRIir[!(dLnLtabPRIir$locname %in% dLnLtabPRIirOUTB$locname),]

#Xenarthra
dLnLtabPRIirOUTX <- dLnLtabPRIir %>% 
  rowwise() %>% mutate(max45 = max(MonoX,Star)) %>%
  arrange(desc(max45)) %>% head(1000) %>% mutate(dif = MonoX-Star)

dLnLtabPRIirNONOUTX <- dLnLtabPRIir[!(dLnLtabPRIir$locname %in% dLnLtabPRIirOUTX$locname),]

#check summed likelihoods for the outlier / non outlier subsets 

#Afrotheria
apply(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirOUTA$locname,c(2,5)], 2, sum, na.rm=T)
apply(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirNONOUTA$locname,c(2,5)], 2, sum, na.rm=T)

#Boreoeutheria
apply(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirOUTB$locname,c(3,5)], 2, sum, na.rm=T)
apply(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirNONOUTB$locname,c(3,5)], 2, sum, na.rm=T)

#Xenarthra
apply(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirOUTX$locname,c(4,5)], 2, sum, na.rm=T)
apply(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirNONOUTX$locname,c(4,5)], 2, sum, na.rm=T)

#AU tests for nonoutliers and outliers separately

#Afrotheria
library(scaleboot)
AUtestPRIiNONOUTA <- relltest(na.omit(as.matrix(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirNONOUTA$locname,c(2,5)])),nb=1000)
AUtestPRIiNONOUTAsummary <- summary(AUtestPRIiNONOUTA)
AUtestPRIiNONOUTAsummary

AUtestPRIiOUTA <- relltest(na.omit(as.matrix(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirOUTA$locname,c(2,5)])),nb=1000)
AUtestPRIiOUTAsummary <- summary(AUtestPRIiOUTA)
AUtestPRIiOUTAsummary

AUtestPRIiNONOUTAsummaryAttr <- attributes(AUtestPRIiNONOUTAsummary)
AUtestPRIiNONOUTApval <- 1-AUtestPRIiNONOUTAsummaryAttr$table$value[which.max(AUtestPRIiNONOUTAsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIiNONOUTApval

AUtestPRIiOUTAsummaryAttr <- attributes(AUtestPRIiOUTAsummary)
AUtestPRIiOUTApval <- 1-AUtestPRIiOUTAsummaryAttr$table$value[which.max(AUtestPRIiOUTAsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIiOUTApval

#Boreoeutheria
AUtestPRIiNONOUTB <- relltest(na.omit(as.matrix(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirNONOUTB$locname,c(3,5)])),nb=1000)
AUtestPRIiNONOUTBsummary <- summary(AUtestPRIiNONOUTB)
AUtestPRIiNONOUTBsummary

AUtestPRIiOUTB <- relltest(na.omit(as.matrix(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirOUTB$locname,c(3,5)])),nb=1000)
AUtestPRIiOUTBsummary <- summary(AUtestPRIiOUTB)
AUtestPRIiOUTBsummary

AUtestPRIiNONOUTBsummaryAttr <- attributes(AUtestPRIiNONOUTBsummary)
AUtestPRIiNONOUTBpval <- 1-AUtestPRIiNONOUTBsummaryAttr$table$value[which.max(AUtestPRIiNONOUTBsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIiNONOUTBpval

AUtestPRIiOUTBsummaryAttr <- attributes(AUtestPRIiOUTBsummary)
AUtestPRIiOUTBpval <- 1-AUtestPRIiOUTBsummaryAttr$table$value[which.max(AUtestPRIiOUTBsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIiOUTBpval

#Xenarthra
AUtestPRIiNONOUTX <- relltest(na.omit(as.matrix(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirNONOUTX$locname,c(4,5)])),nb=1000)
AUtestPRIiNONOUTXsummary <- summary(AUtestPRIiNONOUTX)
AUtestPRIiNONOUTXsummary

AUtestPRIiOUTX <- relltest(na.omit(as.matrix(dLnLtabPRIi_wide[dLnLtabPRIi_wide$V1 %in% dLnLtabPRIirOUTX$locname,c(4,5)])),nb=1000)
AUtestPRIiOUTXsummary <- summary(AUtestPRIiOUTX)
AUtestPRIiOUTXsummary

AUtestPRIiNONOUTXsummaryAttr <- attributes(AUtestPRIiNONOUTXsummary)
AUtestPRIiNONOUTXpval <- 1-AUtestPRIiNONOUTXsummaryAttr$table$value[which.max(AUtestPRIiNONOUTXsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIiNONOUTXpval

AUtestPRIiOUTXsummaryAttr <- attributes(AUtestPRIiOUTXsummary)
AUtestPRIiOUTXpval <- 1-AUtestPRIiOUTXsummaryAttr$table$value[which.max(AUtestPRIiOUTXsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIiOUTXpval

#outlier figure; skipping for now

#separate loci by the topology they favor (have the highest signal for)

#add wiggle room to difference in constraint & star LnL (e.g. >=110%)

#this is just creating the function for the next chunk
get_favoring_lnl_loci <- function(inputdf, signal="dLnL", m){
  inputdf <- select(inputdf, locname, m, Star)
  inputdfb <- apply(inputdf[,2:3], 1, function(x) which(x == max(x)))
  inputdfa <- apply(inputdf[,2:3], 1, function(x) x[which(x == max(x))])
  returndf <- as.data.frame(inputdf$locname[lapply(inputdfa, length) == 1])
  returndf$topology <- as.numeric(unlist(inputdfb[lapply(inputdfb, length) == 1]))
  returndf$dLnL <- as.numeric(unlist(inputdfa[lapply(inputdfa, length) == 1]))
  colnames(returndf) <- c("locname","topology",signal)
  returndf$topology <- factor(returndf$topology, labels = colnames(inputdf)[2:3])
  return(returndf)
}
 
#get groups of loci for each topo (For MonoA v. Star comparison); can't compare across hypotheses, rather only on this binary comparison
dLnLtabPRIirFA <- get_favoring_lnl_loci(dLnLtabPRIir, m = "MonoA")
dLnLtabPRIirFfiltA <- na.omit(dLnLtabPRIirFA)

#get groups of loci for each topo (MonoB vs. Star)
dLnLtabPRIirFB <- get_favoring_lnl_loci(dLnLtabPRIir, m = "MonoB")
dLnLtabPRIirFfiltB <- na.omit(dLnLtabPRIirFB)

#get groups of loci for each topo (MonoX vs. Star)
dLnLtabPRIirFX <- get_favoring_lnl_loci(dLnLtabPRIir, m = "MonoX")
dLnLtabPRIirFfiltX <- na.omit(dLnLtabPRIirFX)

#then parse between outliers & non-outliers

#Afrotheria
dLnLtabPRIirFfiltOUTA <- dLnLtabPRIirFfiltA[dLnLtabPRIirFfiltA$locname %in% dLnLtabPRIirOUTA$locname,]
dLnLtabPRIirFfiltNONOUTA <- dLnLtabPRIirFfiltA[dLnLtabPRIirFfiltA$locname %in% dLnLtabPRIirNONOUTA$locname,]

#Boreoeutheria
dLnLtabPRIirFfiltOUTB <- dLnLtabPRIirFfiltB[dLnLtabPRIirFfiltB$locname %in% dLnLtabPRIirOUTB$locname,]
dLnLtabPRIirFfiltNONOUTB <- dLnLtabPRIirFfiltB[dLnLtabPRIirFfiltB$locname %in% dLnLtabPRIirNONOUTB$locname,]

#Xenarthra
dLnLtabPRIirFfiltOUTX <- dLnLtabPRIirFfiltX[dLnLtabPRIirFfiltX$locname %in% dLnLtabPRIirOUTX$locname,]
dLnLtabPRIirFfiltNONOUTX <- dLnLtabPRIirFfiltX[dLnLtabPRIirFfiltX$locname %in% dLnLtabPRIirNONOUTX$locname,]

#check and test counts; significance of difference in # of loci favoring one topology over the other
library(rstatix)

#Afrotheria
table(dLnLtabPRIirFfiltA$topology)
chisq_all_PRIiA <- chisq_test(as.data.frame(table(dLnLtabPRIirFfiltA$topology))[,2])
chisq_all_PRIiA
Nloci_all_PRIiA <-  chisq_all_PRIiA$p

table(dLnLtabPRIirFfiltOUTA$topology)
chisq_out_PRIiA <- chisq_test(as.data.frame(table(dLnLtabPRIirFfiltOUTA$topology))[,2])
chisq_out_PRIiA
Nloci_out_PRIi <-  chisq_out_PRIiA$p

table(dLnLtabPRIirFfiltNONOUTA$topology)
chisq_nonout_PRIiA <- chisq_test(as.data.frame(table(dLnLtabPRIirFfiltNONOUTA$topology))[,2])
chisq_nonout_PRIiA
Nloci_nonout_PRIiA <-  chisq_nonout_PRIiA$p

#Boreoeutheria
table(dLnLtabPRIirFfiltB$topology)
chisq_all_PRIiB <- chisq_test(as.data.frame(table(dLnLtabPRIirFfiltB$topology))[,2])
chisq_all_PRIiB
Nloci_all_PRIiB <-  chisq_all_PRIiB$p

table(dLnLtabPRIirFfiltOUTB$topology)
chisq_out_PRIiB <- chisq_test(as.data.frame(table(dLnLtabPRIirFfiltOUTB$topology))[,2])
chisq_out_PRIiB
Nloci_out_PRIiB <-  chisq_out_PRIiB$p

table(dLnLtabPRIirFfiltNONOUTB$topology)
chisq_nonout_PRIiB <- chisq_test(as.data.frame(table(dLnLtabPRIirFfiltNONOUTB$topology))[,2])
chisq_nonout_PRIiB
Nloci_nonout_PRIiB <-  chisq_nonout_PRIiB$p

#Xenarthra
table(dLnLtabPRIirFfiltX$topology)
chisq_all_PRIiX <- chisq_test(as.data.frame(table(dLnLtabPRIirFfiltX$topology))[,2])
chisq_all_PRIiX
Nloci_all_PRIiX <-  chisq_all_PRIiX$p

table(dLnLtabPRIirFfiltOUTX$topology)
chisq_out_PRIiX <- chisq_test(as.data.frame(table(dLnLtabPRIirFfiltOUTX$topology))[,2])
chisq_out_PRIiX
Nloci_out_PRIiX <-  chisq_out_PRIiX$p

table(dLnLtabPRIirFfiltNONOUTX$topology)
chisq_nonout_PRIiX <- chisq_test(as.data.frame(table(dLnLtabPRIirFfiltNONOUTX$topology))[,2])
chisq_nonout_PRIiX
Nloci_nonout_PRIiX <-  chisq_nonout_PRIiX$p

#check and test signal strength - Median LnL per constraint/Star, also btwn outliers & non-outliers

#Afrotheria
median(dLnLtabPRIirFfiltA$dLnL[dLnLtabPRIirFfiltA$topology == "MonoA"])
median(dLnLtabPRIirFfiltA$dLnL[dLnLtabPRIirFfiltA$topology == "Star"])

median(dLnLtabPRIirFfiltOUTA$dLnL[dLnLtabPRIirFfiltOUTA$topology == "MonoA"])
median(dLnLtabPRIirFfiltOUTA$dLnL[dLnLtabPRIirFfiltOUTA$topology == "Star"])

median(dLnLtabPRIirFfiltNONOUTA$dLnL[dLnLtabPRIirFfiltNONOUTA$topology == "MonoA"])
median(dLnLtabPRIirFfiltNONOUTA$dLnL[dLnLtabPRIirFfiltNONOUTA$topology == "Star"])

lnLsignal_PRIiA <- run_custom_tests_unpaired(dLnLtabPRIirFfiltA, "dLnL")
lnLsignal_PRIiOUTA <- run_custom_tests_unpaired(dLnLtabPRIirFfiltOUTA,"dLnL")
lnLsignal_PRIiNONOUTA <- run_custom_tests_unpaired(dLnLtabPRIirFfiltNONOUTA, "dLnL")

#Boreoeutheria
median(dLnLtabPRIirFfiltB$dLnL[dLnLtabPRIirFfiltB$topology == "MonoB"])
median(dLnLtabPRIirFfiltB$dLnL[dLnLtabPRIirFfiltB$topology == "Star"])

median(dLnLtabPRIirFfiltOUTB$dLnL[dLnLtabPRIirFfiltOUTB$topology == "MonoB"])
median(dLnLtabPRIirFfiltOUTB$dLnL[dLnLtabPRIirFfiltOUTB$topology == "Star"])

median(dLnLtabPRIirFfiltNONOUTB$dLnL[dLnLtabPRIirFfiltNONOUTB$topology == "MonoB"])
median(dLnLtabPRIirFfiltNONOUTB$dLnL[dLnLtabPRIirFfiltNONOUTB$topology == "Star"])

lnLsignal_PRIiB <- run_custom_tests_unpaired(dLnLtabPRIirFfiltB, "dLnL")
lnLsignal_PRIiOUTB <- run_custom_tests_unpaired(dLnLtabPRIirFfiltOUTB,"dLnL")
lnLsignal_PRIiNONOUTB <- run_custom_tests_unpaired(dLnLtabPRIirFfiltNONOUTB, "dLnL")

#Xenarthra
median(dLnLtabPRIirFfiltX$dLnL[dLnLtabPRIirFfiltX$topology == "MonoX"])
median(dLnLtabPRIirFfiltX$dLnL[dLnLtabPRIirFfiltX$topology == "Star"])

median(dLnLtabPRIirFfiltOUTX$dLnL[dLnLtabPRIirFfiltOUTX$topology == "MonoX"])
median(dLnLtabPRIirFfiltOUTX$dLnL[dLnLtabPRIirFfiltOUTX$topology == "Star"])

median(dLnLtabPRIirFfiltNONOUTX$dLnL[dLnLtabPRIirFfiltNONOUTX$topology == "MonoX"])
median(dLnLtabPRIirFfiltNONOUTX$dLnL[dLnLtabPRIirFfiltNONOUTX$topology == "Star"])

lnLsignal_PRIiX <- run_custom_tests_unpaired(dLnLtabPRIirFfiltX, "dLnL")
lnLsignal_PRIiOUTX <- run_custom_tests_unpaired(dLnLtabPRIirFfiltOUTX,"dLnL")
lnLsignal_PRIiNONOUTX <- run_custom_tests_unpaired(dLnLtabPRIirFfiltNONOUTX, "dLnL")





###Concatenated fit

#Concat data prep

dLnLtabPRIc <- read.csv("combined_iqtree_dLnLs_concat.csv", header = F, stringsAsFactors = F)
dLnLtabPRIc <- distinct(dLnLtabPRIc)
dLnLtabPRIc$V1 <- paste0(dLnLtabPRIc$V1, ".fasta")

amatabfiltPRIc <- amatab[match (dLnLtabPRIc$V1,amatab$Alignment_name ),]

#AU tests

#Afrotheria
AUtestPRIcA <- relltest(na.omit(as.matrix(dLnLtabPRIc[,c(2:5)])),nb=1000)
AUtestPRIcAsummary <- summary(AUtestPRIcA)
AUtestPRIcAsummary

AUtestPRIcAsummaryAttr <- attributes(AUtestPRIcAsummary)
AUtestPRIcApval <- 1-AUtestPRIcAsummaryAttr$table$value[which.max(AUtestPRIcAsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIcApval

#Boreoeutheria
AUtestPRIcB <- relltest(na.omit(as.matrix(dLnLtabPRIc[,c(3:5)])),nb=1000)
AUtestPRIcBsummary <- summary(AUtestPRIcB)
AUtestPRIcBsummary

AUtestPRIcBsummaryAttr <- attributes(AUtestPRIcBsummary)
AUtestPRIcBpval <- 1-AUtestPRIcBsummaryAttr$table$value[which.max(AUtestPRIcBsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIcBpval

#Xenarthra
AUtestPRIcX <- relltest(na.omit(as.matrix(dLnLtabPRIc[,c(4:5)])),nb=1000)
AUtestPRIcXsummary <- summary(AUtestPRIcX)
AUtestPRIcXsummary

AUtestPRIcXsummaryAttr <- attributes(AUtestPRIcXsummary)
AUtestPRIcXpval <- 1-AUtestPRIcXsummaryAttr$table$value[which.max(AUtestPRIcXsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIcXpval

#prep rescaled data
dLnLtabPRIcr <- dLnLtabPRIc
dLnLtabPRIcr[,2] <- as.numeric(dLnLtabPRIcr[,2])
dLnLtabPRIcr[,3] <- as.numeric(dLnLtabPRIcr[,3])
dLnLtabPRIcr[,4] <- as.numeric(dLnLtabPRIcr[,4])
dLnLtabPRIcr[,5] <- as.numeric(dLnLtabPRIcr[,5])
dLnLtabPRIcr <- na.omit(dLnLtabPRIcr)
dLnLtabPRIcr[,2:5] <- (dLnLtabPRIcr[,2:5] - apply(dLnLtabPRIcr[,2:5], 1, min))/amatabfiltPRIc$Alignment_length
colnames(dLnLtabPRIcr) <- c("locname","MonoA","MonoB","MonoX","Star")

#get outlier/nonoutlier subsets

#again, doing it differently than Alex's code

#Afrotheria
dLnLtabPRIcrOUTA <- dLnLtabPRIcr %>% 
  rowwise() %>% mutate(max25 = max(MonoA,Star)) %>%
  arrange(desc(max25)) %>% head(1000) %>% mutate(dif = MonoA-Star)

dLnLtabPRIcrNONOUTA <- dLnLtabPRIcr[!(dLnLtabPRIcr$locname %in% dLnLtabPRIcrOUTA$locname),]

#Boreoeutheria
dLnLtabPRIcrOUTB <- dLnLtabPRIcr %>% 
  rowwise() %>% mutate(max35 = max(MonoB,Star)) %>%
  arrange(desc(max35)) %>% head(1000) %>% mutate(dif = MonoB-Star)

dLnLtabPRIcrNONOUTB <- dLnLtabPRIcr[!(dLnLtabPRIcr$locname %in% dLnLtabPRIcrOUTB$locname),]

#Xenarthra
dLnLtabPRIcrOUTX <- dLnLtabPRIcr %>% 
  rowwise() %>% mutate(max45 = max(MonoX,Star)) %>%
  arrange(desc(max45)) %>% head(1000) %>% mutate(dif = MonoX-Star)

dLnLtabPRIcrNONOUTX <- dLnLtabPRIcr[!(dLnLtabPRIcr$locname %in% dLnLtabPRIcrOUTX$locname),]

#AU Tests on subsets

#Afrotheria
AUtestPRIcNONOUTA <- relltest(na.omit(as.matrix(dLnLtabPRIc[dLnLtabPRIc$V1 %in% dLnLtabPRIcrNONOUTA$locname,c(2:5)])),nb=1000)
AUtestPRIcNONOUTAsummary <- summary(AUtestPRIcNONOUTA)
AUtestPRIcNONOUTAsummary

AUtestPRIcOUTA <- relltest(na.omit(as.matrix(dLnLtabPRIc[dLnLtabPRIc$V1 %in% dLnLtabPRIcrOUTA$locname,c(2:5)])),nb=1000)
AUtestPRIcOUTAsummary <- summary(AUtestPRIcOUTA)
AUtestPRIcOUTAsummary

AUtestPRIcNONOUTAsummaryAttr <- attributes(AUtestPRIcNONOUTAsummary)
AUtestPRIcNONOUTApval <- 1-AUtestPRIcNONOUTAsummaryAttr$table$value[which.max(AUtestPRIcNONOUTAsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIcNONOUTApval

AUtestPRIcOUTAsummaryAttr <- attributes(AUtestPRIcOUTAsummary)
AUtestPRIcOUTApval <- 1-AUtestPRIcOUTAsummaryAttr$table$value[which.max(AUtestPRIcOUTAsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIcOUTApval

#Boreoeutheria
AUtestPRIcNONOUTB <- relltest(na.omit(as.matrix(dLnLtabPRIc[dLnLtabPRIc$V1 %in% dLnLtabPRIcrNONOUTB$locname,c(3:5)])),nb=1000)
AUtestPRIcNONOUTBsummary <- summary(AUtestPRIcNONOUTB)
AUtestPRIcNONOUTBsummary

AUtestPRIcOUTB <- relltest(na.omit(as.matrix(dLnLtabPRIc[dLnLtabPRIc$V1 %in% dLnLtabPRIcrOUTB$locname,c(3:5)])),nb=1000)
AUtestPRIcOUTBsummary <- summary(AUtestPRIcOUTB)
AUtestPRIcOUTBsummary

AUtestPRIcNONOUTBsummaryAttr <- attributes(AUtestPRIcNONOUTBsummary)
AUtestPRIcNONOUTBpval <- 1-AUtestPRIcNONOUTBsummaryAttr$table$value[which.max(AUtestPRIcNONOUTBsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIcNONOUTBpval

AUtestPRIcOUTBsummaryAttr <- attributes(AUtestPRIcOUTBsummary)
AUtestPRIcOUTBpval <- 1-AUtestPRIcOUTBsummaryAttr$table$value[which.max(AUtestPRIcOUTBsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIcOUTBpval

#Xenarthra
AUtestPRIcNONOUTX <- relltest(na.omit(as.matrix(dLnLtabPRIc[dLnLtabPRIc$V1 %in% dLnLtabPRIcrNONOUTX$locname,c(4:5)])),nb=1000)
AUtestPRIcNONOUTXsummary <- summary(AUtestPRIcNONOUTX)
AUtestPRIcNONOUTXsummary

AUtestPRIcOUTX <- relltest(na.omit(as.matrix(dLnLtabPRIc[dLnLtabPRIc$V1 %in% dLnLtabPRIcrOUTX$locname,c(4:5)])),nb=1000)
AUtestPRIcOUTXsummary <- summary(AUtestPRIcOUTX)
AUtestPRIcOUTXsummary

AUtestPRIcNONOUTXsummaryAttr <- attributes(AUtestPRIcNONOUTXsummary)
AUtestPRIcNONOUTXpval <- 1-AUtestPRIcNONOUTXsummaryAttr$table$value[which.max(AUtestPRIcNONOUTXsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIcNONOUTXpval

AUtestPRIcOUTXsummaryAttr <- attributes(AUtestPRIcOUTXsummary)
AUtestPRIcOUTXpval <- 1-AUtestPRIcOUTXsummaryAttr$table$value[which.max(AUtestPRIcOUTXsummaryAttr$table$value[,"sk.3"]),"sk.3"]
AUtestPRIcOUTXpval

#get groups of loci for each topology

#get groups of loci for each topo (For MonoA v. Star comparison); can't compare across hypotheses, rather only on this binary comparison
dLnLtabPRIcrFA <- get_favoring_lnl_loci(dLnLtabPRIcr, m = "MonoA")
dLnLtabPRIcrFfiltA <- na.omit(dLnLtabPRIcrFA)

#get groups of loci for each topo (MonoB vs. Star)
dLnLtabPRIcrFB <- get_favoring_lnl_loci(dLnLtabPRIcr, m = "MonoB")
dLnLtabPRIcrFfiltB <- na.omit(dLnLtabPRIcrFB)

#get groups of loci for each topo (MonoX vs. Star)
dLnLtabPRIcrFX <- get_favoring_lnl_loci(dLnLtabPRIcr, m = "MonoX")
dLnLtabPRIcrFfiltX <- na.omit(dLnLtabPRIcrFX)

#then parse between outliers & non-outliers

#Afrotheria
dLnLtabPRIcrFfiltOUTA <- dLnLtabPRIcrFfiltA[dLnLtabPRIcrFfiltA$locname %in% dLnLtabPRIcrOUTA$locname,]
dLnLtabPRIcrFfiltNONOUTA <- dLnLtabPRIcrFfiltA[dLnLtabPRIcrFfiltA$locname %in% dLnLtabPRIcrNONOUTA$locname,]

#Boreoeutheria
dLnLtabPRIcrFfiltOUTB <- dLnLtabPRIcrFfiltB[dLnLtabPRIcrFfiltB$locname %in% dLnLtabPRIcrOUTB$locname,]
dLnLtabPRIcrFfiltNONOUTB <- dLnLtabPRIcrFfiltB[dLnLtabPRIcrFfiltB$locname %in% dLnLtabPRIcrNONOUTB$locname,]

#Xenarthra
dLnLtabPRIcrFfiltOUTX <- dLnLtabPRIcrFfiltX[dLnLtabPRIcrFfiltX$locname %in% dLnLtabPRIcrOUTX$locname,]
dLnLtabPRIcrFfiltNONOUTX <- dLnLtabPRIcrFfiltX[dLnLtabPRIcrFfiltX$locname %in% dLnLtabPRIcrNONOUTX$locname,]

#check counts and stats

#Afrotheria
table(dLnLtabPRIcrFfiltA$topology)
chisq_all_PRIcA <- chisq_test(as.data.frame(table(dLnLtabPRIcrFfiltA$topology))[,2])
chisq_all_PRIcA
Nloci_all_PRIcA <-  chisq_all_PRIcA$p

table(dLnLtabPRIcrFfiltOUTA$topology)
chisq_out_PRIcA <- chisq_test(as.data.frame(table(dLnLtabPRIcrFfiltOUTA$topology))[,2])
chisq_out_PRIcA
Nloci_out_PRIcA <-  chisq_out_PRIcA$p

table(dLnLtabPRIcrFfiltNONOUTA$topology)
chisq_nonout_PRIcA <- chisq_test(as.data.frame(table(dLnLtabPRIcrFfiltNONOUTA$topology))[,2])
chisq_nonout_PRIcA
Nloci_nonout_PRIcA <-  chisq_nonout_PRIcA$p

#Boreoeutheria
table(dLnLtabPRIcrFfiltB$topology)
chisq_all_PRIcB <- chisq_test(as.data.frame(table(dLnLtabPRIcrFfiltB$topology))[,2])
chisq_all_PRIcB
Nloci_all_PRIcB <-  chisq_all_PRIcB$p

table(dLnLtabPRIcrFfiltOUTB$topology)
chisq_out_PRIcB <- chisq_test(as.data.frame(table(dLnLtabPRIcrFfiltOUTB$topology))[,2])
chisq_out_PRIcB
Nloci_out_PRIcB <-  chisq_out_PRIcB$p

table(dLnLtabPRIcrFfiltNONOUTB$topology)
chisq_nonout_PRIcB <- chisq_test(as.data.frame(table(dLnLtabPRIcrFfiltNONOUTB$topology))[,2])
chisq_nonout_PRIcB
Nloci_nonout_PRIcB <-  chisq_nonout_PRIcB$p

#Xenarthra
table(dLnLtabPRIcrFfiltX$topology)
chisq_all_PRIcX <- chisq_test(as.data.frame(table(dLnLtabPRIcrFfiltX$topology))[,2])
chisq_all_PRIcX
Nloci_all_PRIcX <-  chisq_all_PRIcX$p

table(dLnLtabPRIcrFfiltOUTX$topology)
chisq_out_PRIcX <- chisq_test(as.data.frame(table(dLnLtabPRIcrFfiltOUTX$topology))[,2])
chisq_out_PRIcX
Nloci_out_PRIcX <-  chisq_out_PRIcX$p

table(dLnLtabPRIcrFfiltNONOUTX$topology)
chisq_nonout_PRIcX <- chisq_test(as.data.frame(table(dLnLtabPRIcrFfiltNONOUTX$topology))[,2])
chisq_nonout_PRIcX
Nloci_nonout_PRIcX <-  chisq_nonout_PRIcX$p

#signal strength tests - Median LnL per constraint/Star, also btwn outliers & non-outliers

#Afrotheria
median(dLnLtabPRIcrFfiltA$dLnL[dLnLtabPRIcrFfiltA$topology == "MonoA"])
median(dLnLtabPRIcrFfiltA$dLnL[dLnLtabPRIcrFfiltA$topology == "Star"])

median(dLnLtabPRIcrFfiltOUTA$dLnL[dLnLtabPRIcrFfiltOUTA$topology == "MonoA"])
median(dLnLtabPRIcrFfiltOUTA$dLnL[dLnLtabPRIcrFfiltOUTA$topology == "Star"])

median(dLnLtabPRIcrFfiltNONOUTA$dLnL[dLnLtabPRIcrFfiltNONOUTA$topology == "MonoA"])
median(dLnLtabPRIcrFfiltNONOUTA$dLnL[dLnLtabPRIcrFfiltNONOUTA$topology == "Star"])

lnLsignal_all_PRIcA <- run_custom_tests_unpaired(dLnLtabPRIcrFfiltA, "dLnL")
lnLsignal_out_PRIcA <- run_custom_tests_unpaired(dLnLtabPRIcrFfiltOUTA,"dLnL")
lnLsignal_nonout_PRIcA <- run_custom_tests_unpaired(dLnLtabPRIcrFfiltNONOUTA, "dLnL")

#Boreoeutheria
median(dLnLtabPRIcrFfiltB$dLnL[dLnLtabPRIcrFfiltB$topology == "MonoA"])
median(dLnLtabPRIcrFfiltB$dLnL[dLnLtabPRIcrFfiltB$topology == "Star"])

median(dLnLtabPRIcrFfiltOUTB$dLnL[dLnLtabPRIcrFfiltOUTB$topology == "MonoA"])
median(dLnLtabPRIcrFfiltOUTB$dLnL[dLnLtabPRIcrFfiltOUTB$topology == "Star"])

median(dLnLtabPRIcrFfiltNONOUTB$dLnL[dLnLtabPRIcrFfiltNONOUTB$topology == "MonoA"])
median(dLnLtabPRIcrFfiltNONOUTB$dLnL[dLnLtabPRIcrFfiltNONOUTB$topology == "Star"])

lnLsignal_all_PRIcB <- run_custom_tests_unpaired(dLnLtabPRIcrFfiltB, "dLnL")
lnLsignal_out_PRIcB <- run_custom_tests_unpaired(dLnLtabPRIcrFfiltOUTB,"dLnL")
lnLsignal_nonout_PRIcB <- run_custom_tests_unpaired(dLnLtabPRIcrFfiltNONOUTB, "dLnL")

#Xenarthra
median(dLnLtabPRIcrFfiltX$dLnL[dLnLtabPRIcrFfiltX$topology == "MonoA"])
median(dLnLtabPRIcrFfiltX$dLnL[dLnLtabPRIcrFfiltX$topology == "Star"])

median(dLnLtabPRIcrFfiltOUTX$dLnL[dLnLtabPRIcrFfiltOUTX$topology == "MonoA"])
median(dLnLtabPRIcrFfiltOUTX$dLnL[dLnLtabPRIcrFfiltOUTX$topology == "Star"])

median(dLnLtabPRIcrFfiltNONOUTX$dLnL[dLnLtabPRIcrFfiltNONOUTX$topology == "MonoA"])
median(dLnLtabPRIcrFfiltNONOUTX$dLnL[dLnLtabPRIcrFfiltNONOUTX$topology == "Star"])

lnLsignal_all_PRIcX <- run_custom_tests_unpaired(dLnLtabPRIcrFfiltX, "dLnL")
lnLsignal_out_PRIcX <- run_custom_tests_unpaired(dLnLtabPRIcrFfiltOUTX,"dLnL")
lnLsignal_nonout_PRIcX <- run_custom_tests_unpaired(dLnLtabPRIcrFfiltNONOUTX, "dLnL")




###SCFs 

#prep the data and stats
scf_tabPRI <- read.csv("combined_iqtree_mammals_scf.csv", header = F, stringsAsFactors = F)

#screen_data_filt_scfPRI <- screen_data[match (scf_tabPRI$V1,screen_data$locname ),]
#scf_tabPRI <- scf_tabPRI[screen_data_filt_scfPRI$Rsq>0.25,]
#don't have screened by BLC list

#no SCF values for Star tree (can't make that calculation)
scf_tabPRI$V2 <- as.numeric(scf_tabPRI$V2)
scf_tabPRI$V3 <- as.numeric(scf_tabPRI$V3)
scf_tabPRI$V4 <- as.numeric(scf_tabPRI$V4)
colnames(scf_tabPRI) <- c("locname", "MonoA","MonoB","MonoX")

#Parsing dataframe for each constraint

#Afrotheria
scf_tabPRIA <- select(scf_tabPRI, locname, MonoA) 

#Boreoeutheria
scf_tabPRIB <- select(scf_tabPRI, locname, MonoB)

#Xenarthra
scf_tabPRIX <- select(scf_tabPRI, locname, MonoX)




###Relationship btwn lnL and sCF- For Monophyletic clade constraints

#combining tables; for any particular contig, highest likelihood + topo & value of SCF 

#Individual fit

#Afrotheria
dLnL_VS_sCFtabPRIiA <- merge(dLnLtabPRIirFA, scf_tabPRIA, by="locname")

#Boreoeutheria
dLnL_VS_sCFtabPRIiB <- merge(dLnLtabPRIirFB, scf_tabPRIB, by="locname")

#Xenarthra
dLnL_VS_sCFtabPRIiX <- merge(dLnLtabPRIirFX, scf_tabPRIX, by="locname")

#Concatenated fit

#Afrotheria
dLnL_VS_sCFtabPRIcA <- merge(dLnLtabPRIcrFA, scf_tabPRIA, by="locname")

#Boreoeutheria
dLnL_VS_sCFtabPRIcB <- merge(dLnLtabPRIcrFB, scf_tabPRIB, by="locname")

#Xenarthra
dLnL_VS_sCFtabPRIcX <- merge(dLnLtabPRIcrFX, scf_tabPRIX, by="locname")


###Barplots for LnLs - Within Comparisons

#Individual fit

#Afrotheria
LnLMonoA_Only <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='MonoA', -3]
LnLStar_A <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='Star', -3]

MonoA_Count <- nrow(LnLMonoA_Only)
StarA_Count <- nrow(LnLStar_A)
MonoA_Comp_LnL_Barplot <- barplot(c(MonoA_Count,StarA_Count), names.arg=c('Monophyletic','Star'), ylab='Loci Count', main='Ind. Fit LnL - Afrotheria Comparison', ylim=c(0,80000))

#Boreoeutheria
LnLMonoB_Only <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='MonoB', -3]
LnLStar_B <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='Star', -3]

MonoB_Count <- nrow(LnLMonoB_Only)
StarB_Count <- nrow(LnLStar_B)
MonoB_Comp_LnL_Barplot <- barplot(c(MonoB_Count,StarB_Count), names.arg=c('Monophyletic','Star'),main='Ind. Fit LnL - Boreoeutheria Comparison',ylab='Loci Count',ylim=c(0,80000))

#Xenarthra
LnLMonoX_Only <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='MonoX', -3]
LnLStar_X <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='Star', -3]

MonoX_Count <- nrow(LnLMonoX_Only)
StarX_Count <- nrow(LnLStar_X)
MonoX_Comp_LnL_Barplot <- barplot(c(MonoX_Count, StarX_Count),names.arg=c('Monophyletic','Star'),main='Ind. Fit LnL - Xenarthra Comparison',ylab='Loci Count',ylim=c(0,80000))

#Concatenated fit

#Afrotheria
LnLMonoAc_Only <- dLnLtabPRIcrFA[dLnLtabPRIcrFA$topology=='MonoA', -3]
LnLStar_Ac <- dLnLtabPRIcrFA[dLnLtabPRIcrFA$topology=='Star', -3]

MonoAc_Count <- nrow(LnLMonoAc_Only)
StarAc_Count <- nrow(LnLStar_Ac)
MonoAc_Comp_LnL_Barplot <- barplot(c(MonoAc_Count,StarAc_Count),names.arg=c('Monophyletic','Star'),main='Concat. Fit LnL - Afrotheria Comparison',ylab='Loci Count',ylim=c(0,5e+5))

#Boreoeutheria
LnLMonoBc_Only <- dLnLtabPRIcrFB[dLnLtabPRIcrFB$topology=='MonoB', -3]
LnLStar_Bc <- dLnLtabPRIcrFB[dLnLtabPRIcrFB$topology=='Star', -3]

MonoBc_Count <- nrow(LnLMonoBc_Only)
StarBc_Count <- nrow(LnLStar_Bc)
MonoBc_Comp_LnL_Barplot <- barplot(c(MonoBc_Count,StarBc_Count),names.arg=c('Monophyletic','Star'),main='Concat. Fit LnL - Boreoeutheria Comparison',ylab='Loci Count',ylim=c(0,5e+5))

#Xenarthra
LnLMonoXc_Only <- dLnLtabPRIcrFX[dLnLtabPRIcrFX$topology=='MonoX', -3]
LnLStar_Xc <- dLnLtabPRIcrFX[dLnLtabPRIcrFX$topology=='Star', -3]

MonoXc_Count <- nrow(LnLMonoXc_Only)
StarXc_Count <- nrow(LnLStar_Xc)
MonoXc_Comp_LnL_Barplot <- barplot(c(MonoXc_Count,StarXc_Count),names.arg=c('Monophyletic','Star'),main='Concat. Fit LnL - Xenarthra Comparison',ylab='Loci Count',ylim=c(0,5e+5))

###Barplots for sCFs - Within Comparisons

###Barplot for LnL - Across Comparisons

#Individual fit

#Loci w/ favored topology being respective constraint; also remove dLnL column

LnLMonoA_Only <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='MonoA', -3]

LnLMonoB_Only <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='MonoB', -3]

LnLMonoX_Only <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='MonoX', -3]

#Loci w/ LnL that favor all 3 constraints compared to star

LnLConstraintAB <- merge(LnLMonoA_Only, LnLMonoB_Only, by="locname")
#can only merge two data frames at once
LnLConstraintsFavored <- merge(LnLConstraintAB, LnLMonoX_Only, by="locname")
Constraints_Count <- nrow(LnLConstraintsFavored)

#Loci w/ LnL that have a mix of favored topologies across comparisons (e.g. MonoA, MonoB, Star; or Star, Star, MonoX)
#MonoA, MonoB, Star

LnLMonoA_Only <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='MonoA', -3]

LnLMonoB_Only <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='MonoB', -3]

LnLStar_X <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='Star', -3]

LnLConstraintAB <- merge(LnLMonoA_Only, LnLMonoB_Only, by="locname")
LnL_A_B_Star <- merge(LnLConstraintAB, LnLStar_X, by="locname")
A_B_Star_Count <- nrow(LnL_A_B_Star)

#MonoA, Star, MonoX

LnLMonoA_Only <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='MonoA', -3]

LnLStar_B <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='Star', -3]

LnLMonoX_Only <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='MonoX', -3]

LnLConstraintAX <- merge(LnLMonoA_Only, LnLMonoX_Only, by="locname")
LnL_A_Star_X <-  merge(LnLConstraintAX, LnLStar_B, by="locname")
A_Star_X_Count <- nrow(LnL_A_Star_X)

#Star, MonoB, MonoX

LnLStar_A <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='Star', -3]

LnLMonoB_Only <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='MonoB', -3]

LnLMonoX_Only <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='MonoX', -3]

LnLConstraintBX <- merge(LnLMonoB_Only, LnLMonoX_Only, by="locname")
LnL_Star_B_X <- merge(LnLConstraintBX, LnLStar_A, by="locname")
Star_B_X_Count <- nrow(LnL_Star_B_X)

#Count of 2 Monos, 1 Star - Total
Mono2_Star1_Count <- sum(A_B_Star_Count+A_Star_X_Count+Star_B_X_Count)

#MonoA, Star, Star

LnLMonoA_Only <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='MonoA', -3]

LnLStar_B <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='Star', -3]

LnLStar_X <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='Star', -3]

LnLStarBX <- merge(LnLStar_B, LnLStar_X, by="locname")
LnL_A_Star_Star <- merge(LnLStarBX, LnLMonoA_Only, by="locname")
A_Star_Star_Count <- nrow(LnL_A_Star_Star)

#Star, MonoB, Star

LnLStar_A <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='Star', -3]

LnLMonoB_Only <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='MonoB', -3]

LnLStar_X <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='Star', -3]

LnLStarAX <- merge(LnLStar_A, LnLStar_X, by="locname")
LnL_Star_B_Star <- merge(LnLStarAX, LnLMonoB_Only, by="locname")
Star_B_Star_Count <- nrow(LnL_Star_B_Star)

#Star, Star, MonoX

LnLStar_A <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='Star', -3]

LnLStar_B <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='Star', -3]

LnLMonoX_Only <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='MonoX', -3]

LnLStarAB <- merge(LnLStar_A, LnLStar_B, by="locname")
LnL_Star_Star_X <- merge(LnLStarAX, LnLMonoX_Only, by="locname")
Star_Star_X_Count <- nrow(LnL_Star_Star_X)

#Count of 1 Mono, 2 Stars - Total
Mono1_Star2_Count <- sum(A_Star_Star_Count+Star_B_Star_Count+Star_Star_X_Count)

#Loci w/ LnL that favor Star in all 3 comparisons compared to constraints

LnLStar_A <- dLnLtabPRIirFA[dLnLtabPRIirFA$topology=='Star', -3]

LnLStar_B <- dLnLtabPRIirFB[dLnLtabPRIirFB$topology=='Star', -3]

LnLStar_X <- dLnLtabPRIirFX[dLnLtabPRIirFX$topology=='Star', -3]

LnLStarAB <- merge(LnLStar_A, LnLStar_B, by="locname")
LnLStarsFavored<- merge(LnLStarAB, LnLStar_X, by="locname")
Stars_Count <- nrow(LnLStarsFavored)

#Vector of category counts & names
Category_Counts <- c(Constraints_Count, Mono2_Star1_Count, Mono1_Star2_Count, Stars_Count)
Category_Names <- c(3, 2, 1, 0)

IndividualFitBarPlot_LnLs <- barplot(Category_Counts, names.arg=Category_Names, main='Individual Fit LnL-Favored Topologies', ylim=c(0,65000), xlab='Monophyletic Constraints Favored', ylab='Loci Count')

#Concatenated fit

#Loci w/ LnL that favor all 3 constraints compared to star

#Loci w/ LnL that have a mix of favored topologies across comparisons (e.g. MonoA, MonoB, Star; or Star, Star, MonoX)

#Loci w/ LnL that favor Star in all 3 comparisons compared to constraints



###Barplot for sCFs - Across Comparisons

#Loci w/ sCF >50 across all 3 comparisons

SCF_MonoA <- scf_tabPRIA[scf_tabPRIA$MonoA>50, ]
SCF_MonoB <- scf_tabPRIB[scf_tabPRIB$MonoB>50, ]
SCF_MonoX <- scf_tabPRIX[scf_tabPRIX$MonoX>50, ]




#Loci w/ mix of sCF values (e.g. MonoA sCF=56, MonoB=34, MonoX=67)

#Loci w/ sCF <50  across all 3 comparisons

