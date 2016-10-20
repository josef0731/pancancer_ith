#MSI CNV analysis

file = "C:/Users/joseph/Desktop/Level_3/msi.txt"
samples <- read.table(file, sep = "\t", header = FALSE, row.names = NULL, 
                  stringsAsFactors = FALSE, quote = "", comment.char = "")

msi = samples[,2]

msi = paste("C:/Users/joseph/Desktop/Level_3/msi/",msi,sep="")

mean_chr7 = rep(1, length(msi))
mean_chr8 = rep(1, length(msi))

for (k in 1:length(msi)){
  
  seg <- read.table(file = msi[k], sep = "\t", header = TRUE, row.names = NULL, 
                    stringsAsFactors = FALSE, quote = "", comment.char = "")
  
  copy_number = 2^(seg[,6] + 1)
  
  seg = cbind(seg, copy_number)
  seg_chr7 = seg[seg[,"Chromosome"] == 7,]
  seg_chr8 = seg[seg[,"Chromosome"] == 8,]
  seg_chr1 = seg[seg[,"Chromosome"] == 1,]
  seg_chr20 = seg[seg[,"Chromosome"] == 20,]
  
  length_chr7 = 159138663
  length_chr8 = 145138636
  length_chr20 = 63025520
  length_chr1 = 249250621
  
  prop_chr7 = (seg_chr7[,"End"] - seg_chr7[,"Start"] + 1) / length_chr7
  prop_chr8 = (seg_chr8[,"End"] - seg_chr8[,"Start"] + 1) / length_chr8
  
  weighted_cn7 = seg_chr7$copy_number * prop_chr7
  weighted_cn8 = seg_chr8$copy_number * prop_chr8
  mean_chr7[k] = sum(weighted_cn7)
  mean_chr8[k] = sum(weighted_cn8)
}

file = "C:/Users/joseph/Desktop/Level_3/TCGA_MSI_CN_chr7_8.txt"

samples = cbind(samples, mean_chr7, mean_chr8)
write.table(samples, file,quote = FALSE, sep="\t", row.names = FALSE)