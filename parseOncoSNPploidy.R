args <- commandArgs(trailingOnly = TRUE)

summaryfile = args[1]

masterfile = read.table(summaryfile, header = TRUE, sep = "\t", skip = 1
                        stringsAsFactors = FALSE, quote = "", comment.char = "",row.names = NULL)

ploidy = masterfile[1,"Copy Number (Average)"]

file = args[2]

case = unlist(strsplit(summaryfile, split = "."))[1]

output = cbind(case, ploidy)

append = TRUE
coln = FALSE
if (!file.exists(file)){
  file.create(file)
  append = FALSE
  coln = TRUE
}

write.table(x = output, file = file, append=append,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names=coln)
