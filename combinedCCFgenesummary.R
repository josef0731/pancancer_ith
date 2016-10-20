combinedCCF = read.table("combinedCCF.txt", header = TRUE, sep = "\t",
                         quote = "", comment.char = "",row.names = NULL)

genes = function(x = 0){
  #MSI drivers according to our Nat Genet paper 2014 (n = 37)
  natgenet14 = c("CASP5","AIM2","OR52N5","IGSF9B","RNF43","SLC16A4","SYNJ2","PRR11","CLOCK","ACVR2A",
                 "KIAA1009","SMAP1","CDC7","SEC63","LMAN1","UVRAG","RBM27","RGS12","HMMR","ZMYM4",
                 "RBM6","KRAS","RBM45","UPF3A","TTK","SLC35F5","MBD6","TFAM","ARID1A","SEC31A",
                 "TEAD2","EBPL","COBLL1","DYNC1I2","TGFBR2","RPL22","TMBIM4")
  
  #Cancer genome census (n = 572)
  cgc = read.csv("C:/Users/joseph/Downloads/Census_allFri Dec  4 04-14-46 2015.csv", header = TRUE,
                 stringsAsFactors = FALSE, comment.char = "", row.names=NULL)
  genes = cgc[,1]
  
  #ccf estimates on BOTH our MSI drivers and CGC census genes (Total: 604)
  genes = union(natgenet14, genes)
  return(genes)
}

genes = genes()

combinedCCF = combinedCCF[which(combinedCCF[,"gene"] %in% genes),]

subtype = seq(1, nrow(combinedCCF))

subtypes <- function(i){
  return(combined[combined[,"case"] == combinedCCF[i,"case"],"subtype"])
}

subtype = sapply(subtype, subtypes)

combinedCCF = cbind(combinedCCF, subtype)

cCCFlist = split.data.frame(combinedCCF, combinedCCF[,"subtype"], drop = TRUE)

for (j in 2:6){
  library(plyr)
  genesummary = ddply(cCCFlist[[j]], .variables = "gene", summarise,
                      clonal = sum(pr_0.85 >= 0.5),
                      subclonal = sum(pr_0.85 < 0.5))
  
  genesummary = genesummary[which(genesummary[,"gene"] %in% genes),]
  
  total = rep("A", nrow(genesummary))
  total = genesummary[,"clonal"] + genesummary[,"subclonal"]
  genesummary = cbind(genesummary, total)
  genesummary = genesummary[order(genesummary$total, decreasing = TRUE),]
  
  
  fisher.pval = rep(1, nrow(genesummary))
  for (i in 1:nrow(genesummary)){
    gene = c(genesummary[i, "clonal"], genesummary[i, "subclonal"])
    rival = c(sum(genesummary[, "clonal"]) - genesummary[i, "clonal"], 
              sum(genesummary[, "subclonal"]) - genesummary[i, "subclonal"])
    fisher.pval[i] = fisher.test(rbind(as.numeric(gene), as.numeric(rival)))$p.value
  }
  
  genesummary = cbind(genesummary, fisher.pval)
  
  source("C:/Users/joseph/Desktop/qval_calc.R")
  qval = qvalue(genesummary[,"fisher.pval"], alpha = 0.05)$qvalue
  
  genesummary = cbind(genesummary, qval)
  
  genesummary = genesummary[order(genesummary$qval, decreasing = FALSE),]
  
  write.table(genesummary, paste("combinedCCFgenesummary", cCCFlist[[j]][1,"subtype"], ".txt", sep = ""),
              quote = FALSE, row.names = FALSE,
              col.names = TRUE, sep = "\t")
}

library(plyr)
genesummary = ddply(combinedCCF, .variables = "gene", summarise,
                    clonal = sum(pr_0.85 >= 0.5),
                    subclonal = sum(pr_0.85 < 0.5))

genesummary = genesummary[which(genesummary[,"gene"] %in% genes),]

total = rep("A", nrow(genesummary))
total = genesummary[,"clonal"] + genesummary[,"subclonal"]
genesummary = cbind(genesummary, total)
genesummary = genesummary[order(genesummary$total, decreasing = TRUE),]


fisher.pval = rep(1, nrow(genesummary))
for (i in 1:nrow(genesummary)){
  gene = c(genesummary[i, "clonal"], genesummary[i, "subclonal"])
  rival = c(sum(genesummary[, "clonal"]) - genesummary[i, "clonal"], 
            sum(genesummary[, "subclonal"]) - genesummary[i, "subclonal"])
  fisher.pval[i] = fisher.test(rbind(as.numeric(gene), as.numeric(rival)))$p.value
}

genesummary = cbind(genesummary, fisher.pval)

source("C:/Users/joseph/Desktop/qval_calc.R")
qval = qvalue(genesummary[,"fisher.pval"], alpha = 0.05)$qvalue

genesummary = cbind(genesummary, qval)

genesummary = genesummary[order(genesummary$qval, decreasing = FALSE),]

write.table(genesummary, "combinedCCFgenesummary.txt",
            quote = FALSE, row.names = FALSE,
            col.names = TRUE, sep = "\t")
