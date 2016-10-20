files <- list.files(pattern = "\\.maf$")

output = cbind("case", "nMut", "tcontent")

for (i in 1:length(files)){
  maf = read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = NULL)
  nMut = nrow(maf)
  maf = maf[maf[,"Mutation.Type"] == "snp",]
  coverage = (maf[,"t_alt_count"] + maf[,"t_ref_count"])
  tcontent = 2 * maf[,"t_alt_count"] / coverage
  tcontent = median(tcontent, na.rm = TRUE)
  if (tcontent > 1) tcontent = 1
  case = files[i]
  print(case)
  app = cbind(case, nMut, tcontent)
  output = rbind(output, app)
}

outputloc = "coverageTcontent.txt"

write.table(output,file = outputloc, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
