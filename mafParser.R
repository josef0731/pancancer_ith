files <- list.files(pattern = "\\.maf$")

for (i in 1:length(files)){
  maf = read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = NULL)
  maf = maf[maf[,"Variant_Classification"] != "RNA",]
  maf[!duplicated(maf), ]
  print(nrow(maf))
  tumor_depth = maf[, "t_alt_count"] + maf[, "t_ref_count"]
  maf = cbind(maf, tumor_depth)
  maf[, "tumor_depth"] = as.numeric(maf[,"tumor_depth"])
  maf = maf[maf[, "tumor_depth"] >= 10, ]
  maf = maf[maf[, "t_alt_count"] >= 3, ]
  newfilename = paste(unlist(strsplit(files[i], split = ".maf"))[1], "parsed.maf", sep = "")
  write.table(maf,file = newfilename, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}
