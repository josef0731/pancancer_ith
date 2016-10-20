masterfile = read.table("C:/Users/joseph/Documents/UCECcnvMasterFile.txt", header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, quote = "", comment.char = "",row.names = NULL)

samples = as.factor(masterfile[,"sample"])

for (i in 1:length(levels(samples))){
  SampleID = as.character(levels(samples)[i])
  print(SampleID)
  TumourFile = masterfile[masterfile[,"sample"] == SampleID & masterfile[,"type"] == "Tumor", "filename"]
  NormalFile = masterfile[masterfile[,"sample"] == SampleID & masterfile[,"type"] == "BloodN", "filename"]
  output = cbind(SampleID, TumourFile, NormalFile)
  if(ncol(output) != 3){
    NormalFile = masterfile[masterfile[,"sample"] == SampleID & masterfile[,"type"] == "SolidN", "filename"]
    output = cbind(SampleID, TumourFile, NormalFile)
  }
  if(ncol(output) == 3){
    if (i != 1) final = rbind(final, output)
    else final = output
  }
}
filename = "C:/Users/joseph/Documents/UCEC_oncosnp_masterfile.txt"
write.table(x = final, file = filename, append = FALSE, quote = FALSE,
            col.names = TRUE, row.names = FALSE, sep = "\t")

oncosnpbatch = read.table("C:/Users/joseph/Documents/UCEC_oncosnp_masterfile.txt", header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, quote = "", comment.char = "",row.names = NULL)

for (i in 1:30){
  if(nrow(oncosnpbatch) < 18*i) small = nrow(oncosnpbatch) else small = 18*i
  output = oncosnpbatch[c(seq((i-1)*18 + 1, small)),]
  filename = paste("C:/Users/joseph/Documents/UCEC_oncosnp_masterfile", i, ".txt", sep = "")
  write.table(x = output, file = filename, append = FALSE, quote = FALSE,
              col.names = TRUE, row.names = FALSE, sep = "\t")
  
}
