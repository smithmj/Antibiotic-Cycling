setwd("~/Antibiotic Cycling") 
source(file="cycling.R")

for(i in 1:length(genotypes)){
  print(genotypes[i])
  initial_state <- rep(0, 16)
  initial_state[i] <- 1
  adaptive_list <- adaptive_cycling(initial_state, 1000)
  class(adaptive_list)
  adaptive_seq <- adaptive_list[[1]]
  seq_growth_rates <- as.numeric(adaptive_list[[2]])
  final_state <- adaptive_list[[3]]
}


                                 