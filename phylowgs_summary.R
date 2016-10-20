args <- commandArgs(trailingOnly = TRUE)

tree_summary = read.csv("top_tree_stat.txt",
                        header=TRUE, sep = "\t")

setting = args[1] #setting: if only protein-altering mutations enter "WES". if whole-genome mutations enter "WGS". 
sample = args[2]

library(plyr)
sum_stat <- ddply(tree_summary,.variables = c("rank","clone."),summarise, 
                  mean_freq=mean(freq), 
                  median_ssm = median(X.ssm), 
                  median_cnv = median(X.cnv), parent = as.integer(median(parent)))

#filter away ssm < 1% of all SSMs and cnv == 0

n_total_mutation = sum(sum_stat[sum_stat[,"rank"] == 0,]$median_ssm)
n_cluster = nrow(sum_stat[sum_stat[,"rank"] == 0,])
sum_stat = sum_stat[sum_stat$median_cnv > 0 | sum_stat$median_ssm >= sum(sum_stat[sum_stat[,"rank"] == 0,]$median_ssm) * 0.01,]
n_new =  sum(sum_stat[sum_stat[,"rank"] == 0,]$median_ssm)
n_mutation_filtered = n_total_mutation - n_new

#inspect sum_stat
#sum_sum_stat <- ddply(sum_stat, .variables = c("clone."), summarise,
#                      mean_freq = mean(mean_freq),
#                      median_ssm = median(median_ssm),
#                      median_cnv = median(median_cnv))

# FOR NOW TAKE ONLY THE MOST FREQUENT TREE. TO IMPROVE NEEDS PARSING MERGING OF THE TREES.
sum_sum_stat = sum_stat[sum_stat[,"rank"] == 0,]
sum_sum_stat = sum_sum_stat[,2:ncol(sum_sum_stat)]
colnames(sum_sum_stat)[1] = "clone"

n_cluster_filtered = n_cluster - nrow(sum_sum_stat)

sum_ssm = sum(sum_sum_stat[,3])
clones = sum_sum_stat[,1]

#parse cnv
cnv_coordinate =  read.csv("../cnv_data.txt",
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
length = cnv_coordinate[,4] - cnv_coordinate[,3]
total = 3098825702
cnv_prop = length / total
cnv_coordinate = cbind(cnv_coordinate[,1:4], cnv_prop)

clone =  read.csv("../top_trees/clone_member.txt",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cnv = rep("X", nrow(clone))
retain = rep(NULL, 0)
for (i in 1:nrow(clone)){
  cnv[i] =  unlist(strsplit(clone[i,2], split = ""))[1]
  if (cnv[i] == "c") retain = append(retain, i)
}

cnv_clone =  clone[retain,]

clone_id = rep(0, nrow(cnv_coordinate))
for (i in 1:length(clone_id)){
  for (j in 1:nrow(cnv_clone)){
    if (cnv_clone[j, 2] == cnv_coordinate[i, 1]) clone_id[i] = cnv_clone[j, 1]
  }
}
cnv_coordinate = cbind(cnv_coordinate, clone_id)


#check proper assignment of cnv. IF no, the script stops and require manual assignment.
continue = TRUE
proportion = rep(0, nrow(sum_sum_stat))
for (k in 1:length(proportion)){
  if (sum(cnv_clone[,1] == clones[k]) != sum_sum_stat[sum_sum_stat[,"clone"] == clones[k],"median_cnv"]){
    continue = FALSE
  }
}

continue = TRUE

#functions for plotting. Load first.
require(shape)
col2 <- shadepalette("limegreen", "black", n = 50)    #no mut

height <- function (clone, mode){
  if (mode == "CNV") index = "proportion"
  else index = "ssm_prop"
  current = clone
  current_height = 0
  while (current != 0){
    cur_parent = parent[current]
    if (cur_parent == 0) break
    current_height = current_height + sum_sum_stat[sum_sum_stat[,1] == cur_parent,index]
    current = cur_parent
  }
  return(current_height - 0.65)
}

center <- function (clone){
  if (clone == 0){
    return(0)
  } 
  if (sum(parent == parent[clone]) < 2) {
    return (center(parent[clone]))
  }
  else{
    current = clone
    current_center = 0
    cur_parent = parent[current]
    if (cur_parent != 0){
      if (sum_sum_stat[min(sum_sum_stat[sum_sum_stat[,"parent"] == parent[clone],1]), 1] == clone){
        current_center = center(cur_parent) + ((sum_sum_stat[sum_sum_stat[,1] == cur_parent,2] / pi) ^ 0.5) * 0.7
      }
      else{
        current_center = center(cur_parent) - ((sum_sum_stat[sum_sum_stat[,1] == cur_parent,2] / pi) ^ 0.5) * 0.7
      }
    }
    if (current_center == 0){
      if (clone != 1){
        current_center = - (sum_sum_stat[clone,2] / pi) ^ (0.5) - 0.05
      }
      else current_center = (sum_sum_stat[clone,2] / pi) ^ (0.5) + 0.05
    } 
  }
  return(current_center)
}

color <- function(mut){
  col = col2
  return (col)
}


draw <- function(mode, setting){
  if(mode == "CNV") {
    emptyplot(c(-0.5, 0.5), c(-0.7, 0.7), cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, main = paste(sample,"CNV", sep = ""), frame.plot = TRUE)
    legend("topleft", legend = paste("Proportion of \n genome with\n CNV:\n ", 
                                   round(prop_CNV,digits = 4), sep = ""), cex = 1.5, bty="n")
  }
  else {
    if (setting == "WES") column = 2
    else column = 3
    emptyplot(c(-0.5, 0.5), c(-0.7, 0.7),  cex.lab=2.5, cex.axis=2.5, cex.main=2.5, cex.sub=2.5, main = paste(sample, "SNV_indel", sep = ""), frame.plot=TRUE)
    record = read.table("../ssm_data.txt",
                        header=TRUE, sep = "\t", stringsAsFactors = FALSE)  #NEED TO CHANGE!
    legend("topleft", legend = paste("#SNV/indel:\n", 
                                   nrow(record), sep = ""), bty="n", cex = 1.5)
  }
  for (clone in 1:nrow(sum_sum_stat)){
    if (mode == "CNV"){
      length = sum_sum_stat[clones[clone], "proportion"]
      col = col2 
    }
    else{
      length = sum_sum_stat[clones[clone], "ssm_prop"]
      col = color(mutation[clone])
    }
    c = center(clone)
    h = height(clone, mode)
    filledcylinder(rx = 0.02, ry = (sum_sum_stat[clones[clone],2] / pi)^.5, angle = 90, 
                   len = length, col = c(col, rev(col)), 
                   mid = c(c, h + length / 2), topcol = col[25], lcol = "black")
  }
}

if (continue == TRUE){
  for (k in 1:length(proportion)){
    proportion[k] = sum(cnv_coordinate[cnv_coordinate[,"clone_id"] == k,5])
  }
  prop_CNV = sum(proportion)
  proportion = proportion / prop_CNV
  
  mutation = rep(0, nrow(sum_sum_stat))
  
#  for (k in 1:length(mutation)){
#    mutation[k] = parse_mutation(k)
#  }
  
  ssm_prop = sum_sum_stat$median_ssm / sum_ssm
  sum_sum_stat = cbind(sum_sum_stat, ssm_prop, proportion, mutation)
  
  case = rep(sample, nrow(sum_sum_stat))
  sum_sum_stat = cbind(case, sum_sum_stat)
  sum_sum_stat
  
  file = "../../HKU66phylowgs.txt"
  
  append = TRUE
  coln = FALSE
  if (!file.exists(file)){
    file.create(file)
    append = FALSE
    coln = TRUE
  }
  
  write.table(x = sum_sum_stat, file = file, append=append,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names=coln)

  file = "../../HKU66phylowgsFiltering.txt"
  
  append = TRUE
  coln = FALSE
  if (!file.exists(file)){
    file.create(file)
    append = FALSE
    coln = TRUE
  }
  
  filter = data.frame(sample=sample, n_cluster_filtered=n_cluster_filtered, n_mutation_filtered = n_mutation_filtered)

  write.table(x = filter, file = file, append=append,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names=coln)

  
  #plotting
  require(shape)
  
  sum_sum_stat = sum_sum_stat[,2:ncol(sum_sum_stat)]
  mutation = sum_sum_stat[,"mutation"]
  parent = sum_sum_stat[,"parent"]
  type = c("SNVindel", "CNV")

  for (j in 1:length(type)){
    file = paste("../../plots/",sample,"_", type[j], ".pdf", sep = "")
    pdf(file, width=6.45, height=3.95)
    draw(type[j], setting)
    dev.off()
  }
}