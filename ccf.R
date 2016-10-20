# cancer cell fraction calculation


 ##########  ''''''''' REMEMBER DELETING THE ORIGINAL FILE BEFORE RE-GENERATING IT.  '''''''''  ##########

#cases = a vector of Sample ID (with pfgXXX for HKU, TCGA-XX-XXXX for TCGA data)
#Tcontent by median mutation - used for HKU data
#Absolute Tcontent estimate: Rank 1 estimate - used for TCGA data

setup <- function(dataset){
  if (dataset == "TCGA"){
    cases = read.table("coverageTcontent.txt", header = TRUE,
                       stringsAsFactors = FALSE, comment.char = "", row.names=NULL)
    case = cases[,1]
    purity = cases[,3]
    vcfs <- paste("maf/", case, "T.maf",sep="")
    cnvs <- paste("oncosnp/cnv/", case, "_oncosnp.txt", sep = "")
  }
  output = list(case, purity, vcfs, cnvs)
  names(output) = c("case", "purity", 'vcfs', "cnvs")
  return(output)
}

parsing = function(dataset, id, vcfs, cnvs){
  vcf <- read.table(vcfs[id], header = TRUE, sep = "\t",stringsAsFactors = FALSE,
                    quote = "", comment.char = "",row.names = NULL)
  
  cnv <- read.table(cnvs[id], header = TRUE, sep = "\t",stringsAsFactors = FALSE,
                    quote = "", comment.char = "",row.names = NULL)
    
  maf = vcf[vcf[,"Mutation_Type"] == "snp",]
  coverage = (maf[,"t_alt_count"] + maf[,"t_ref_count"])
  tcontent = 2 * maf[,"t_alt_count"] / coverage
  purity = median(tcontent, na.rm = TRUE)

  if (dataset == "TCGA"){
    gene_col = "Hugo_Symbol"
    vcf = vcf[which(vcf[,"Chromosome"] %in% c(seq(1,22), "X", "Y")),]
    for (j in 1:ncol(vcf)){
      if(colnames(vcf)[j] == "Start_Position"){
        colnames(vcf)[j] = "Start_position"
      }
      if(colnames(vcf)[j] == "End_Position"){
        colnames(vcf)[j] = "End_position"
      }

      if(colnames(vcf)[j] == "Reference_Allele"){
        colnames(vcf)[j] = "Germline_allele"
      }
      if(colnames(vcf)[j] == "Tumor_Seq_Allele2"){
        colnames(vcf)[j] = "Mutant_allele"
      }
      if(colnames(vcf)[j] == "cds.coordinate"){
        colnames(vcf)[j] = "cDNA_Change"
      }
      if(colnames(vcf)[j] == "Protein.coordinate"){
        colnames(vcf)[j] = "Protein_Change"
      }

    }
  }

  output = list(vcf, cnv, purity)
  names(output) = c("vcf", "cnv", "purity")
  return(output)
}

ccf_calc = function(dataset, case, vcf, cnv, purity){
  #q = tumor Copy number
  #ccf = cancer cell fraction
  #ci_low
  #ci_high
  
  if (dataset == "TCGA"){
    basic = vcf[,c("Chromosome", "Start_position", "End_position", "Germline_allele", "Mutant_allele", 
                    "Mutation_Type", "Sift.prediction", "Sift.score", "Polyphen.prediction", "Polyphen.score", 
		      "Protein.domains", "MA.Uniprot", "MA.variant", "MA.Func..Impact", "MA.FI.score")]
    case = rep(case, nrow(vcf))
    basic = cbind(case, basic)
    start_col = "Start_position"
    end_col = "End_position"
    alt_col = "t_alt_count"
    cov_col = c("t_ref_count", "t_alt_count")
    gene_col = vcf[,"Hugo_Symbol"]
    DNA_description = vcf[,"cDNA_Change"]
    protein_description = vcf[,"Protein_Change"]
  }

  coverage = rep(NULL, nrow(vcf))
  for (i in 1:nrow(vcf)) coverage[i] = sum(vcf[i,cov_col])
  q = rep(NULL, nrow(vcf))
  ccf = rep(NULL, nrow(vcf))
  var = rep(NULL, nrow(vcf))
  ci_low = rep(0.01, nrow(vcf))
  ci_high = rep(0.01, nrow(vcf))
  error_low = rep(NULL, nrow(vcf))
  error_high = rep(NULL, nrow(vcf))
  pr_0.85 = rep(NULL, nrow(vcf))
  
  for (i in 1:nrow(vcf)){
    # obtain tumor Copy number
    q[i] = 2
    for (j in 1:nrow(cnv)){
      if (vcf[i,"Chromosome"] == "X") q[i] = 2
      else if (vcf[i,"Chromosome"] == "Y") q[i] = 2
      else if (as.integer(vcf[i,"Chromosome"]) == cnv[j,"chrom"] && 
               vcf[i,start_col] >= cnv[j, "start"] && 
               vcf[i,end_col] <= cnv[j, "end"]){
        q[i] = cnv[j,"major_cn"] + cnv[j, "minor_cn"]
      }
    }

    c = seq(from = 0.01, to = 1, by = 0.01)
    f_c = purity * c / (2*(1 - purity) + q[i]*purity)
    p_c = dbinom(x = vcf[i,alt_col], size = coverage[i], prob = f_c)

    if (is.na(as.numeric(sum(p_c))) == FALSE & sum(p_c) != 0) {
      norm_pc = p_c / sum(p_c)
      cum_pc = rep(NULL, 100)
      for (m in 1:100) cum_pc[m] = sum(norm_pc[1:m])
      
      ccf[i] = sum(norm_pc * c)
      pr_0.85[i] = sum(norm_pc[86:100])
      
      var[i] = sum(norm_pc * (c - ccf[i])^2)
      
      ci_high[i] = ccf[i] + qnorm(0.975) * (var[i] / coverage[i]) ^ (0.5)
      ci_low[i] = ccf[i] + qnorm(0.025) * (var[i] / coverage[i]) ^ (0.5)
      
      if (ci_high[i] > 1) ci_high[i] = 1
      error_high[i] = ci_high[i] - ccf[i]
      error_low[i] = ccf[i] - ci_low[i]
    }
    else {
      ccf[i] = 0
      pr_0.85[i] = 0
      var[i]
      ci_high[i] = 0
      ci_low[i] = 0
      error_high[i] = 0
      error_low[i] = 0
    }
    if (i == 200) print("Finished variant #200")
    if (i == 500) print("Finished variant #500")
    if (i == 1000) print("Finished variant #1000")

  }  
  print("ccf calculation completed")
  table = data.frame(basic, coverage, gene_col, DNA_description, protein_description,
                     ccf, ci_low, ci_high, error_low, error_high, pr_0.85)
  return(table)
}

#real execution here
dataset = setup("TCGA") #[or TCGA]
for (k in 1:length(dataset$case)){
#call the above functions ...
  print(dataset$case[k])
  sample = dataset$case[k]
  dt = parsing("TCGA", k, dataset$vcfs, dataset$cnvs)
  table = ccf_calc("TCGA", sample, dt$vcf, dt$cnv, dt$purity)
  
  dest = paste("ccf/", sample, "_ccf.txt", sep = "")
  write.table(x = table, file = dest, append=FALSE,
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  
  coverage = table[,"coverage"]
  ccf = table[, "ccf"]
  clonalProp = length(ccf[ccf >= 0.85]) / length(ccf)
  ccf_reciprocal = 1 / ccf
  ccf_reciprocal = sort(ccf_reciprocal)
  t1 = 0.1
  t1v = ccf_reciprocal[t1 * length(ccf_reciprocal)]
  t2 = 0.95
  t2v = ccf_reciprocal[t2 * length(ccf_reciprocal)]
  
  if (length(ccf) > 5000){
    subsample = sample(ccf, 5000)
  } else {
    subsample = ccf
  }
  
  if (length(ccf) >= 3) {
    normality.test = shapiro.test(subsample)
    normality.pval = round(normality.test$p.value, 6)
    normality.stat = round(normality.test$statistic, 6)
  } else {
    normality.pval = "NA"
    normality.stat = "NA"
  }
  slope = round((t2 - t1) / (t2v - t1v), 6)
  euclid_distance = round((ccf ^ 2 + coverage ^ 2) ^ (0.5), 6)
  dispersion = round(var(euclid_distance) / mean(euclid_distance), 6)
  print("stats completed")
  output = cbind(sample, slope, normality.stat, normality.pval, dispersion, clonalProp)
  
  append = TRUE
  file = "ccf/HKU100_stat.txt"
  if (!file.exists(file)){
    file.create(file)
    append = FALSE
  }

  write.table(x = output, file = file, append=append,
              quote = FALSE, sep = "\t", row.names = FALSE,col.names = !append)

}