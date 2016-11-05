setwd("~/Antibiotic Cycling") 
source(file="cycling.R")
source(file="adaptive_cycling/expanded_adaptive_cycling.R")

predict_cycle <- function(seq, initial_state, gr=growth_rates, ab_list = antibiotics){
  seq_growth_rates <- c()
  state <- initial_state
  for(i in 1:500){
    for(ab in seq){
      state <- state %*% cpm(ab, gr=gr)
      growth_rate <- expected_gr(state, ab, gr=gr)
      seq_growth_rates <- c(seq_growth_rates, growth_rate)
    }
  }
  # We want to find the mean of the  growth rates over the
  # last repetition of the cycle.
  l <- length(seq_growth_rates)
  m <- length(seq)
  start_index <- l-m+1
  stop_index <- l
  seq_avg_gr <- mean(seq_growth_rates[start_index:stop_index])
  # now we will see if any monotherapy of abs have lower equilibrium growthrates
  # than the sequence 
  better_abs <- c() 
  for(ab in ab_list){
    M_ab <- cpm(ab, gr=gr)
    M_ab_eq <- eq_mat(M_ab)
    ab_gr <- expected_gr(initial_state, ab, gr=gr)
    if(ab_gr < seq_avg_gr){
      better_abs <- c(better_abs, ab)
    }
  }
  print(paste("seq_avg_gr:", seq_avg_gr))
  print("better monotherapies:")
  print(better_abs)
  return()
}