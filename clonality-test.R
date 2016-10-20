genes = function(x = 0){
  #MSI drivers according to our Nat Genet paper 2014 (n = 37)
  natgenet14 = c("CASP5","AIM2","OR52N5","IGSF9B","RNF43","SLC16A4","SYNJ2","PRR11","CLOCK","ACVR2A",
                 "KIAA1009","SMAP1","CDC7","SEC63","LMAN1","UVRAG","RBM27","RGS12","HMMR","ZMYM4",
                 "RBM6","KRAS","RBM45","UPF3A","TTK","SLC35F5","MBD6","TFAM","ARID1A","SEC31A",
                 "TEAD2","EBPL","COBLL1","DYNC1I2","TGFBR2","RPL22","TMBIM4")
  
  #Cancer genome census (n = 572)
  cgc = read.csv("Census_allFri Dec  4 04-14-46 2015.csv", header = TRUE,
                 stringsAsFactors = FALSE, comment.char = "", row.names=NULL)
  genes = cgc[,1]
  
  #ccf estimates on BOTH our MSI drivers and CGC census genes (Total: 604)
  genes = union(natgenet14, genes)
  return(genes)
}

gene = genes()

file = "HKU34_ccf_SNVCNV.txt"
#CCF estimates of a total of 570 genes selected as described previously

hkuccf = read.table(file = file, header=TRUE, sep="\t", ,stringsAsFactors = FALSE,
                 quote = "", comment.char = "",row.names = NULL)
hkuccf = hkuccf[which(ccf[,"gene"] %in% gene),]

file = "TCGA289_ccf_SNVCNV.txt"
tcgaccf = read.table(file = file, header=TRUE, sep="\t", ,stringsAsFactors = FALSE,
                 quote = "", comment.char = "",row.names = NULL)
tcgaccf = tcgaccf[which(tcgaccf[,"gene"] %in% gene),]

subtype = function(case, mastertable){
  return(as.character(mastertable[mastertable[,"sample"] == case, "subtype"]))
}

clonal_status = function(type, ccf = 0, pr_ccf = 0){
  if (type == "CNV"){
    if (ccf >= 0.85) return("CLONAL") else return("SUBCLONAL")
  } else {
    if (pr_ccf >= 0.5) return("CLONAL") else return ("SUBCLONAL")
  }
}
#molgp = rep(2, nrow(hkuccf))
#for (i in 1:length(molgp)) molgp[i] = subtype(hkuccf[i,"case"], hku)

ccf = rbind(hkuccf, tcgaccf)

clonal = rep(2, nrow(ccf))
for (i in 1:length(clonal)) clonal[i] = clonal_status(ccf[i,"type"], ccf[i, "ccf"], ccf[i, "pr_0.85"])

ccf = cbind(ccf, clonal)
colnames(tcgaccf)[ncol(ccf)] = c("clonal_status")

intestinal = ccf[ccf[,"subtype"] == "Intestinal",]
diffuse = ccf[ccf[,"subtype"] == "Diffuse",]
ebv = ccf[ccf[,"subtype"] == "EBV",]
mixed = ccf[ccf[,"subtype"] == "Mixed",]
msi = ccf[ccf[,"subtype"] == "MSI",]

listed= list(intestinal, diffuse, ebv, mixed, msi)

for (q in 1:length(listed)){
  data = listed[[q]]
  pair = NULL
  
  library(plyr)
  sum_stat <- ddply(msi,.variables = c("case","gene"), summarize,
                    clonal_count = sum(clonal_status == "CLONAL"),
                    subclonal_count = sum(clonal_status == "SUBCLONAL"))
  
  dichot <- ddply(sum_stat,.variables = c("case","gene"), summarize,
                  clonal = if (clonal_count > 0) 1 else 0,
                  subclonal = if (subclonal_count > 0) 1 else 0)
  
  dichot.split = split.data.frame(dichot, dichot$case)
  
  clonal_pair = NULL
  subclonal_pair = NULL
  order = NULL
  reverse = NULL   
  
  pair = data.frame(clonal_pair, subclonal_pair, order, reverse)
  
  for (i in 1:length(dichot.split)){
    df = dichot.split[[i]]
    clonal = df[df[,"clonal"] == 1,"gene"]
    subclonal = df[df[,"subclonal"] == 1,"gene"]
    for (j in 1:length(clonal)){
      for (k in 1:length(subclonal)){
        appended = FALSE
        clonal_pair = clonal[j]
        subclonal_pair = subclonal[k]
        if (nrow(pair) != 0){
          for (m in 1:nrow(pair)){
            if (pair[m,"clonal_pair"] == clonal_pair & pair[m,"subclonal_pair"] == subclonal_pair){
              pair[m,"order"] = as.integer(pair[m,"order"]) + 1
              appended = TRUE
            } else if (pair[m,"subclonal_pair"] == clonal_pair & pair[m,"clonal_pair"] == subclonal_pair){
              pair[m,"reverse"] = as.integer(pair[m,"reverse"]) + 1
              appended = TRUE
            } 
          }        
        }
        if (appended == FALSE){
          order = 1
          reverse = 0
          entry = cbind(clonal_pair, subclonal_pair, order, reverse)
          pair = rbind(pair, entry)
          pair$order = as.integer(pair$order)
          pair$reverse = as.integer(pair$reverse)
          appended = TRUE
        }
      }
    }
  }
  
  subtypes = c("intestinal", "diffuse", "ebv", "mixed", "msi")
  
  write.table(pair, file = paste("pair_", subtypes[q], ".txt", sep = ""), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
}


#total_clonal = sum(dichot$clonal == 1)
#total_subclonal = sum(dichot$subclonal == 1)

#genecount <- ddply(dichot, .variables = "gene", summarize,
#                   clonal = sum(clonal == 1),#number of cases subclonal
#                   subclonal = sum(subclonal == 1))# number of cases subclonal

#casecount <- ddply(dichot, .variables = "case", summarize,
#                  clonal = sum(clonal == 1),# number of genes clonal
#                   subclonal = sum(subclonal == 1))# number of genes subclonal
#rownames(casecount) = casecount[,1]


#total_by_case_clonal = rep(11, nrow(dichot))
#total_by_case_subclonal = rep(11, nrow(dichot))

#for (i in 1:nrow(dichot)){
#  total_by_case_clonal[i] = casecount[dichot$case[i],2]
#  total_by_case_subclonal[i] = casecount[dichot$case[i],3]
#}

#dichot = cbind(dichot, total_by_case_clonal, total_by_case_subclonal)

#genecasecount <- ddply(dichot, .variables = c("case","gene"), summarize,
#                   out_degree = if (clonal == 1) total_by_case_subclonal - subclonal else 0,
#                   in_degree = if (subclonal == 1) total_by_case_clonal - clonal else 0)

#degreecount <- ddply(genecasecount, .variables = "gene", summarize,
#                     total_in_degree = sum(in_degree),
#                     total_out_degree = sum(out_degree),
#                     binom_pval = if(total_out_degree + total_in_degree != 0) binom.test(total_out_degree, total_out_degree + total_in_degree)$p.value else 0)

#degreecount = cbind (degreecount, genecount$clonal, genecount$subclonal)
#colnames(degreecount)[5:6] = c("no_of_cases_clonal", "no_of_cases_subclonal")
#degreecount = degreecount[!(degreecount[,"total_out_degree"] == 0 & degreecount[,"total_in_degree"] == 0),]
#source("C:/Users/joseph/Desktop/qval_calc.R")
#qval = qvalue(degreecount$binom_pval, 0.05)
#degreecount = cbind(degreecount, qval)

file = "C:/Users/joseph/Documents/Combined_EBV_ccf.txt"
#write.table(degreecount, file = file, quote = FALSE, sep = "\t", row.names = FALSE)

df = read.table(file, sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
total_cases = df$no_of_cases_clonal + df$no_of_cases_subclonal
df = cbind(df, total_cases)
df = df[order(decreasing = TRUE, df$total_cases),]
df = df[df$total_cases >= 5,]

pval = rep("A", nrow(df))

for (i in 1:nrow(df)){
  gene = df[i, "gene"]
  nongene = paste("non",gene, sep = "")
  gene_clonal = df[i,"no_of_cases_clonal"]
  gene_subclonal = df[i,"no_of_cases_subclonal"]
  nongene_clonal = sum(df[df$gene != gene,"no_of_cases_clonal"])
  nongene_subclonal = sum(df[df$gene != gene,"no_of_cases_subclonal"])
  fisher_data = data.frame(cbind(c(gene_clonal, gene_subclonal),c(nongene_clonal, nongene_subclonal)), 
                           row.names = c(gene, nongene), stringsAsFactors = FALSE)
  colnames(fisher_data) = c("clonal", "subclonal")
  pval[i] = fisher.test(fisher_data, alternative = "two.sided")$p.value
}

pval = as.numeric(pval)

source("C:/Users/joseph/Desktop/qval_calc.R")
qval = qvalue(pval, alpha = 0.05)$qvalue

clonal_status = cbind(df[,c("gene", "no_of_cases_clonal", "no_of_cases_subclonal", "total_cases")], pval, qval)

write.table(clonal_status, "clonal_status_EBV.txt", quote = FALSE, sep = "\t", row.names = FALSE)

file = "ccf.txt"
ccf = read.table(file, sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
ccf = split(ccf, ccf$subtype)

for (i in 1:length(ccf)){
  subtype = ccf[[i]][1,"subtype"]
  hist(ccf[[i]][,"ccf"], main = subtype)
}