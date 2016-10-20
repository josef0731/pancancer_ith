sim_NR <- function(sp, sd, mu, Tp, Td, k, binding, g){
  
  # Event 1 = Birth (time = 1 / B(d,p))
  # Event 2 = Death (time = 1 / D(d,p), set = 0 initially)
  # Event 3 = Gain Driver (= mu * Td)
  # Event 4 = Gain Passenger (= mu * Tp)
  # Event 5 = Lose Passenger through immunity 
  # (= (1 - (1 - binding) ^ total_mut) * total_mut / 3 * 10^-9), i.e.
  #       % of genome that is mutated and at least one mut is binding)

  #Figs 1B, 2B, 2C
  
  generation = 0
  generation_array = generation
  
  d = 0 #driver
  p = 0 #passenger
  N = k #population size; initial set at k
  N_array = N
  D_array = d
  P_array = p
  cell_list = rep(list(c(d, p)), N)

  birth_array = rep(1, N)
  mean_birth = mean(birth_array)

  #alphas
  gain_driver = mu * Td #rate to gain driver
  gain_passenger = mu * Tp #rate to gain passenger

  while(TRUE){
    to_be_deleted = NULL
    new_cell_list = NULL
    for (i in 1:length(cell_list)){
      d = cell_list[[i]][1]
      p = cell_list[[i]][2]
      birth = (1 + sd)^d / (1 + sp)^p
      death = N / k # death rate
      if (N >= 10^6) death = log(1 + N / k)
      rate = c(rpois(1, birth), rpois(1, death))
      choice = which(rate %in% max(rate))

      if (2 %in% choice) {
        to_be_deleted = c(to_be_deleted, i)
        N = N - 1
        N_array = c(N_array, N)
        generation_array = c(generation_array, generation)
      } else {
        for (j in 1:2){
          driver_increment = rpois(1, gain_driver)
          passenger_increment = rpois(1, gain_passenger)
          new_cell_list = c(new_cell_list, list(c(d + driver_increment, p + passenger_increment)))
          birth_array = c(birth_array, (1 + sd)^(d + driver_increment) / (1 + sp)^(p + passenger_increment))
        }
        N = N + 2
        N_array = c(N_array, N)
        generation_array = c(generation_array, generation)
      } 
    }
    if (!is.null(to_be_deleted)){
      cell_list = cell_list[-to_be_deleted]
    }
    if (!is.null(new_cell_list)){
      cell_list = c(cell_list, new_cell_list)
    }
    
    print(length(cell_list))
    
    mean_birth = mean(birth_array)
    
    generation = generation + 1 / mean_birth
    
    if (N >= 2 * k | N <= 0 | generation >= g) {
      if (N >= 2 * k) print("SUCCESS")
      break
    }
  }
  return(list(N = N_array, cell = cell_list, generation = generation_array))
}

pdf("simulation_default.pdf", width = 11.69, height = 8.27)
for (i in 1:5){
  if (i != 1) {
    par(new = TRUE)
    xlab = ""
    ylab = ""
    axes = FALSE
    main = ""
  } else {
    xlab = "Time"
    ylab = "Population"
    axes = TRUE
    main = "Default"
  }
  sim = sim_NR(0.001, 0.1, 10^-8, 5*10^6, 700, 1000, 0, 15000)
  plot(y = sim$N, x = sim$generation, xlim = c(0, 15000), ylim = c(0, 2000), 
       xlab = xlab, ylab = ylab, axes = axes, main = main)
}
dev.off()