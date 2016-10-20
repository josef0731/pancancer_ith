masterfile = read.table("../newcaselist.txt", stringsAsFactors = FALSE, header = FALSE, sep = "\t")
cases = masterfile[, 1]
files <- paste(cases, "parsedTRIMMED.maf", sep = "")

output = cbind("case", "n_Mut")

for (i in 1:length(files)){
  maf = read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = NULL)
  case = files[i]
  print(case)
  n_Mut = nrow(maf)
  app = cbind(case, n_Mut)
  output = rbind(output, app)
}

outputloc = "n_MutNEW.txt"

write.table(output,file = outputloc, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
