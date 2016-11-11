setwd("~/Antibiotic Cycling") 
source(file="cycling.R")

#for(i in 1:length(genotypes)){
#  print(genotypes[i])
#  initial_state <- rep(0, 16)
#  initial_state[i] <- 1
#  adaptive_list <- adaptive_cycling(initial_state, 1000)
#  class(adaptive_list)
#  adaptive_seq <- adaptive_list[[1]]
#  seq_growth_rates <- as.numeric(adaptive_list[[2]])
#  final_state <- adaptive_list[[3]]
#}

final_state_mat <- matrix(0, nrow=1000, ncol=16)
avg_growth_rates <- c()
random_initial_states <- read.table("adaptive_cycling/tables/random_initial_states_1000.txt", sep=",")
for(i in 1:nrow(random_initial_states)){
  initial_state <- as.matrix(random_initial_states[i,])
  print(initial_state)
  adaptive_list <- adaptive_cycling(initial_state, 1000)
  adaptive_seq <- adaptive_list[[1]]
  seq_growth_rates <- as.numeric(adaptive_list[[2]])
  l <- length(seq_growth_rates)
  avg_gr <- seq_growth_rates[l-100:l]
  avg_growth_rates <- c(avg_growth_rates, avg_gr)
  final_state <- adaptive_list[[3]]
  print(final_state)
  final_state_mat[i,] <- final_state
}

                                 