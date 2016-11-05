setwd("~/Antibiotic Cycling") 
source(file="cycling.R")

random_states <- as.matrix(read.table("adaptive_cycling/tables/random_initial_states_1000.txt", sep=",", nrows=1000))

sequences <- c()
avg_growth_rates <- c()
max_growth_rates <- c()
seq_len <- 300
final_state_matrix <- matrix(0, nrow=seq_len, ncol=16)
for(i in 1:nrow(random_states)){
  print(i)
  adaptive_list <- adaptive_cycling(random_states[1,], seq_len)
  adaptive_seq <- adaptive_list[[1]]
  sequences <- c(sequences, adaptive_seq)
  seq_growth_rates <- as.numeric(adaptive_list[[2]])
  avg_growth_rates <- c(avg_growth_rates, mean(seq_growth_rates[seq_len-100:seq_len]))
  max_growth_rates <- c(max_growth_rates, max(seq_growth_rates[seq_len-100:seq_len]))
  final_state <- adaptive_list[[3]]
  print(final_state)
  final_state_matrix[i,] <- final_state
}

sequences_table <- as.table(matrix(sequences, nrow=seq_len))
print("mean growth rates:")
print(mean(avg_growth_rates))
print("max growth rates:")
print(mean(max_growth_rates))
avg_final_state <- colMeans(final_state_matrix)
print(avg_final_state)
final_state_sd <- apply(final_state_matrix, 2, sd)
print(final_state_sd)

#for each state, need to save the equilibrium growth rate over the last 100,
#the max growth rate over the last 100, and the sequence so that you
#can find the onset step manually later