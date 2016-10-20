pyclone_id <- read.table(file = "pyclone_id.txt", 
                  sep = "\t", header = TRUE, row.names = NULL, 
                  stringsAsFactors = FALSE, quote = "", comment.char = "")

cellfreq <- read.table(file = "cellular.frequencies.tsv", 
                         sep = "\t", header = TRUE, row.names = NULL, 
                         stringsAsFactors = FALSE, quote = "", comment.char = "")
cellfreq = cellfreq[1001:10000,]

sample = unlist(strsplit(getwd(), split = "/"))[8]

mean_cellfreq = rep(1, ncol(cellfreq))

for (i in 1:length(mean_cellfreq)){
  mean_cellfreq[i] = mean(cellfreq[,i])
}

id = colnames(cellfreq)

cellfreq = data.frame(id, mean_cellfreq)

pyclone_id = pyclone_id[order(pyclone_id$pyclone_id),]

library(plyr)
occurences <- count(pyclone_id,vars = "pyclone_id")

occurences = occurences[order(occurences$freq, decreasing = TRUE),]
n_cluster_filtered = nrow(occurences)

occurences = occurences[occurences$freq > sum(occurences$freq) * 0.01,]
n_cluster_filtered = n_cluster_filtered - nrow(occurences)
n_mutations_filtered = nrow(pyclone_id) - sum(occurences$freq)

k = nrow(occurences)
sum = sum(occurences$freq)
prop = occurences$freq / sum
max_prop = max(prop)
simpson = 1 / sum(prop * prop)

logprop = rep("s", k)
for (i in 1:k) if (prop[i] != 0) logprop[i] = log(prop[i]) else prop[i] = 0
logprop = as.numeric(logprop)

shannon = - sum(logprop * prop)

output = cbind(sample, k, max_prop, shannon, simpson, n_cluster_filtered, n_mutations_filtered)

file = "../../HKU66_pyclone_summary.txt"

append = TRUE
coln = FALSE
if (!file.exists(file)){
  file.create(file)
  append = FALSE
  coln = TRUE
}

write.table(x = output, file = file, append=append,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names=coln)

cellfreq = cellfreq[order(cellfreq$id),]
pyclone_id = pyclone_id[order(pyclone_id$mutation),]
pyclone_id = cbind(pyclone_id, cellfreq$mean_cellfreq)
write.table(pyclone_id, "pyclone_id_new.txt", sep = "\t",
            quote = FALSE, row.names = FALSE,col.names = TRUE)