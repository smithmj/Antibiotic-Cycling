setwd("~/Antibiotic Cycling") 
source(file="cycling.R")

random_initial_states <- read.table("adaptive_cycling/tables/random_initial_states_1000.txt", sep=",")
num_cycles <- nrow(random_initial_states)

### adaptive cycling ###
cycle_length <- 1000 # must be greater than 100
final_state_mat <- matrix(0, nrow=num_cycles, ncol=16)
avg_growth_rates <- c()
adaptive_seq_mat <- matrix(0, nrow=num_cycles, ncol=cycle_length)
for(i in 1:nrow(random_initial_states)){
  print(i)
  initial_state <- as.matrix(random_initial_states[i,])
  adaptive_list <- adaptive_cycling(initial_state, cycle_length)
  adaptive_seq <- adaptive_list[[1]]
  adaptive_seq_mat[i,] <- adaptive_seq
  seq_growth_rates <- adaptive_list[[2]]
  avg_gr <- seq_growth_rates[cycle_length-100:cycle_length]
  avg_growth_rates <- c(avg_growth_rates, avg_gr)
  final_state <- adaptive_list[[3]]
  print(final_state)
  final_state_mat[i,] <- final_state
}

# find the average equilibrium distribution across all the initial states
colMeans(final_state_mat)
# find the standard deviation of the equilbrium distribution across all the initial states
apply(final_state_mat, 2, sd)

# find the average growth rate over all the initial states
mean(avg_growth_rates)
# find the standard deviation of the average growth rates over all the initial states
sd(avg_growth_rates)

### AMC monotherapy ###
AMC_final_state_mat <- matrix(0, nrow=num_cycles, ncol=16)
for(i in 1:nrow(random_initial_states)){
  initial_state <- as.matrix(random_initial_states[i,])
  AMC_eq_state <- initial_state %*% eq_mat(cpm("AMC"))
  AMC_final_state_mat[i,] <- AMC_eq_state
}
# find the average equilibrium distribution across all the initial states
colMeans(AMC_final_state_mat)
# find the standard deviation of the equilbrium distribution across all the initial states
apply(AMC_final_state_mat, 2, sd)



