require(mcclust)
labels = read.table("labels.tsv", header = TRUE, sep="\t")
labels = labels[1001:10000,]
labels = labels + 1
labels = as.matrix(labels)

psm <- comp.psm(labels)
# posterior similarity matrix
# optimize criteria based on PSM

mpear <- maxpear(psm)
# Relabelling
mutation = colnames(labels)
pyclone_id = mpear$cl
output = cbind (mutation, pyclone_id)

write.table(output, "pyclone_id.txt", quote=FALSE,
            sep = "\t", row.names = FALSE)

psm = data.frame(psm)
psm = round(psm, digits = 6)
colnames(psm) = mutation
rownames(psm) = mutation

#reorder psm according to pyclone id
id = output
id = id[order(id[,"pyclone_id"]),]
mutation = id[,"mutation"]
psm = psm[match(mutation,rownames(psm)),]
psm = psm[,match(mutation,colnames(psm))]

write.table(psm, "pyclone_psm.txt", quote=FALSE, sep = "\t")