msi = c("057","103","127","129","143","212")

clone_members <- paste(msi,"_clone_member.txt",sep="")
vcfs <- paste("pfg",msi,"T",sep="")

for (k in 1:1){
  vcf <- read.table(vcfs[k], header = TRUE, sep = "\t",stringsAsFactors = FALSE,
                    quote = "", comment.char = "",row.names = NULL)
  clone_member <- read.table(clone_members[k], header = TRUE, sep = "\t",stringsAsFactors = FALSE,
                              quote = "", comment.char = "",row.names = NULL)
  
  chr = rep(1, nrow(clone_member))
  pos = rep(1, nrow(clone_member))
  
  for (i in 1:nrow(clone_member)){
    chr[i] = as.numeric(unlist(strsplit(clone_member$id, split = "_")[i])[1])
    pos[i] = as.numeric(unlist(strsplit(clone_member$id, split = "_")[i])[2])
  }
  
  #strip away "chr" from chromosome numbers
  
  strip = function(x){
    x = unlist(strsplit(x, split = "r"))[2]
    #    chr = try(as.integer(chr),silent = TRUE)
  }  
  
  vcf[,2] = sapply(vcf[,2], function(x) unlist(strsplit(x,split = "r"))[2])
  vcf = vcf[vcf[,2] != "X",]
  vcf = vcf[vcf[,2] != "Y",]
  
  clone_member = cbind(clone_member, chr, pos)
  clone_member = clone_member[order(clone_member$chr,clone_member$pos),]
  mutation_type = rep("N",nrow(clone_member))
  vcf = vcf[order(as.numeric(vcf[,2]),vcf[,3]),]

  for (i in 1:nrow(clone_member)){
    
    if (i == 1) start = 1
    
    for (j in start:nrow(vcf)){
      if (clone_member[i,"chr"] == as.numeric(vcf[j,2]) && 
          any(clone_member[i,"pos"] == seq(vcf[j,3] - 2, vcf[j,3] + 2, by = 1))) {
        mutation_type[i] = vcf[j,7]
        start = j
        break
      }
    }
  }
  
  clone_member = cbind (clone_member, mutation_type)
  
  file = paste(msi[k],"_clone_member_with_mutation_type.txt",sep="")
  
  write.table(x = clone_member, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
}

