ccf = read.table("C:/Users/joseph/Documents/msi_ccf_MSIdrivers.txt",header = TRUE,sep = "\t")
#ccf = ccf[ccf[,"description_protein"] != "-",]
ccf = ccf[order(ccf$gene,ccf$ccf),] #order data frame by group then by value

gene = ccf$gene
class = rep("UNCLASSIFIED", nrow(ccf))
colour = rep(176,nrow(ccf))
error_colour = rep(229,nrow(ccf))

for (i in 1:nrow(ccf)){
#  if (ccf$ci_high[i] == 1){
  if (ccf$ci_high[i] == 1 & ccf$pr[i] >= 0.75){
    class[i] = "CLONAL"
    colour[i] = 552
    error_colour[i] = 101
  } 
#  if (ccf$ci_high[i] != 1){
  if (ccf$ci_high[i] != 1 & ccf$pr[i] <= 0.25){
    class[i] = "SUBCLONAL"
    colour[i] = 132
    error_colour[i] = 62
  } 
}

cancer_cell_fraction = ccf$ccf

df = data.frame(gene, cancer_cell_fraction, class, colour, error_colour)
df = df[order(df$gene,df$cancer_cell_fraction),] #order data frame by group then by value

#transform group information to surrogate x value for x axis plotting

x = rep(NULL, 0)

summary = summary(df$gene)
summary = summary[summary != 0]

for (i in 1:length(attr(summary,which = "names"))){
  x[length(x) + 1 : summary[i]] = seq(i - 0.75, i - 0.25, length.out = summary[i])
}

labels = attr(summary, which="names")

par(las = 2) #labels perpendicular to axis
palette(colors())

plot(x = x, y = df$cancer_cell_fraction,type="n",bty="n",
     ylim = c(0,1), ylab = "Cancer cell fraction", xlab = "", xaxt = "n", pch = 20, col = df$colour)

arrows(x, ccf$ci_low, x, ccf$ci_high, length = 0, col = df$error_colour)

points(x = x, y = df$cancer_cell_fraction,type="p",ann=FALSE,bty="n",
     ylim = c(0,1), xaxt = "n", pch = 20, col = df$colour)

Axis(x = seq(0,length(summary),1), 
     at = seq(0.5,length(summary) - 0.5,1),
     side=1, labels=labels)

legend("bottomright", bty = "n", pch = 20, legend=c("Clonal","Unclassified","Subclonal"), col=c(552,176,132),
      pt.cex = 1.2, y.intersp =  0.5)
