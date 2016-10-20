arguments <- commandArgs(trailingOnly = TRUE)
case = arguments[1]
gene_col_name = arguments[2]
ccf_file = arguments[3]
driver_list = arguments[4]
case_list = arguments[5]

col <- function(ref, t1, t2){
  if (ref == "-" | t1 == "-" | t2 == "-") color = "red"
  else color = 1
  return(color)
}

plot_ccf <- function(case, gene_col_name, ccf_file, driver_list, case_list){
  
  file = paste(case, "ccfplot_only_nonsynonymous.pdf", sep = "_")
  pdf(file, height = 8.27, width = 11.69)
  
  ccf = read.table(ccf_file, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  
  drivers = read.table(driver_list, stringsAsFactors = FALSE, sep = "\t", header = FALSE)
  drivers = drivers[,1]
  print(drivers)

  cases = read.table(case_list, stringsAsFactors = FALSE, sep = "\t", header = FALSE)
  cases = cases[,1]

#  ccf = ccf[which(ccf$case %in% cases),]

#  plot(y = ccf[,"ccf"], x = seq(1, length(ccf[,"ccf"])), type = "n", main = case, ylab = "ccf",
#       xlab = "mutation #", ylim = c(0, 1), cex.main = 3, cex.axis = 2.6, cex.lab = 2.6)
  table = ccf[,c("ccf", "type", gene_col_name, "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")]
  table = table[table[,"type"] != "Silent",]
  table = table[order(table$ccf, decreasing = TRUE),]
  table[,"ccf"] = as.numeric(table[,"ccf"])

  plot(y = table[,"ccf"], x = seq(1, length(table[,"ccf"])), type = "n", main = case, ylab = "ccf",
       xlab = "mutation #", ylim = c(0, 1), cex.main = 3, cex.axis = 2.6, cex.lab = 2.6)

  for(p in 1:nrow(table)){
    color = col(table[p, "Reference_Allele"], table[p, "Tumor_Seq_Allele1"], table[p, "Tumor_Seq_Allele2"])
    points(x = p, y = table[p, "ccf"], pch = 19, col = color)
    if (as.character(table[p, gene_col_name]) %in% drivers & table[p, "type"] != "Silent"){
      segments(p, y0 = table[p, "ccf"] - 0.15, y1 = table[p, "ccf"] + 0.03, lwd = 2, col = color)
      text(x = p, y = table[p, "ccf"] - 0.15, labels = table[p, gene_col_name], col = color,
           cex = 2, pos = 2, srt = 90)
    }
  }
  
  dev.off()
  
  return("All done")  
}

plot_ccf(case = case, gene_col_name = gene_col_name, ccf_file = ccf_file, 
         driver_list = driver_list, case_list = case_list)


