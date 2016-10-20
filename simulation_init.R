library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

parallelFunc <- function(x){
  source("simulation_nextReactionORIG.R")
  sim_NR(sp = 0.001, sd = 0.1, mu = 10^(-6), Tp = 5 * 10^6, 
         Td = 700, k = 1000, binding = 0.8, g = 15000, clonality = TRUE,
         clonality_plot = paste("clonality_", x, ".pdf", sep = ""), 
         N_plot = paste("Nplot_", x, ".pdf", sep = ""), 
         alternative = TRUE)
  return(paste("done", x, sep = " "))
}

parLapply(cl, 1:5, parallelFunc)