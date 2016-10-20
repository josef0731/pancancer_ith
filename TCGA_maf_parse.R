BI_maf <- read.table(file = "HKU100_proteinAltering_maf.txt", 
                     sep = "\t", header = TRUE, row.names = NULL, 
                     stringsAsFactors = FALSE, quote = "", comment.char = "#")

cn = colnames(BI_maf)

BI_maf_split = split.data.frame(x = BI_maf, f = BI_maf[,"Tumor_Sample_Barcode"])

for (i in 1:length(BI_maf_split)){
  sampleID = attributes(BI_maf_split[i])$names
  print(sampleID)
  filename = paste("maf/", sampleID,".maf", sep = "")
  df = as.data.frame(BI_maf_split[i])
  colnames(df) = cn
  write.table(x = df, file = filename, append = FALSE, quote = FALSE,
              col.names = TRUE, row.names = FALSE, sep = "\t")
}