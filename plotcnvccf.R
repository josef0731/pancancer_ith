a = rnorm(n = 50, mean = 0.7, sd = 0.03)
b = rnorm(n = 30, mean = 0.5, sd = 0.05)
c = rnorm(n = 20, mean = 0.3, sd = 0.1)
hist(a, border = "white", xlim = c(-1,0), ylim = c(0, 15), 
     xaxt = 'n', yaxt = "n", xlab = "", ylab = "", main = "")
lines(density(c, adjust = 2), col="green", lwd=2)
axis(side = 1, at = seq(-1,0,by = 0.1), labels = FALSE, tck = -0.01)
hist(a, border = "white", xlim = c(0,1), ylim = c(0, 15), 
     xaxt = 'n', yaxt = "n", xlab = "", ylab = "", main = "")
lines(density(a, adjust = 2), col="blue", lwd=2)
lines(density(b, adjust = 2), col="red", lwd=2)
axis(side = 1, at = seq(0,1,by = 0.1), labels = FALSE, tck = -0.01)

df = df[order(df$group, df$val),]

abline(v = c(-0.85, 0.85), col = "brown", lwd = 2)
abline(v = 0, col = "black")

subtype = function(case, mastertable){
  return(mastertable[mastertable[,"sample"] == case, "subtype"])
}
color = function(subtype){
  if(subtype == "MSI") return("red")
  else if(subtype == "Intestinal") return("blue")
  else if(subtype == "Diffuse") return("gold2")
  else if(subtype == "Mixed") return("darkolivegreen3")
  else if(subtype == "EBV") return("magenta3")
  else return(NULL)
}

subtypes = c("MSI", "Intestinal", "Diffuse", "Mixed", "EBV")

for (chr in 1:22){
  df = summary[summary[,"Chromosome"] == chr,]
  msi
  intestinal
  diffuse
  mixed
  ebv
  meanlist = list(msi, intestinal, diffuse, mixed, ebv)
  #copy loss first
  hist(a, border = "white", xlim = c(-1,0), ylim = c(0, 15), 
       xaxt = 'n', yaxt = "n", xlab = "", ylab = "", main = "")
  axis(side = 1, at = seq(-1,0,by = 0.1), labels = FALSE, tck = -0.01)
  for (i in 1:5){
    if(mean(list[i]) < 0){
      list[i] = -list[i] 
      lines(density(list[i], adjust = 2), col=color(subtypes[i]), lwd=2)
    } 
  }
  abline(v = c(-0.85, 0.85), col = "brown", lwd = 2)
  
  #copy gain now
  hist(a, border = "white", xlim = c(0,1), ylim = c(0, 15), 
       xaxt = 'n', yaxt = "n", xlab = "", ylab = "", main = "")
  axis(side = 1, at = seq(0,1,by = 0.1), labels = FALSE, tck = -0.01)
  for (i in 1:5){
    if(mean(list[i]) > 0){
      lines(density(list[i], adjust = 2), col=color(subtypes[i]), lwd=2)
    } 
  }
  abline(v = c(-0.85, 0.85), col = "brown", lwd = 2)
  
}