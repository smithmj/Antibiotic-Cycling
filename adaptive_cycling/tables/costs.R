setwd("~/Antibiotic Cycling")
source(file="cycling.R")

# given a vector, this function finds the convergence index n, defined as the minimum value n where vec[m] = vec[n] 
# for all m >= n
find_convergence_step <- function(vec){
  for(n in 1:length(vec)){
    val <- vec[n]
    index_found <- TRUE
    for(m in n:length(vec)){
      if(vec[m] != vec[n]){
        index_found <- FALSE
        break
      }
    }
    if(index_found){
      return(n)
    }
  }
}

find_costs <- function(initial_state, monotherapy_ab, adaptive_seq_length = 1000){
  ## find longterm cost
  # find monotherapy equilibrium growth rate
  state <- initial_state 
  ab_grs <- c()
  for(i in 1:adaptive_seq_length){
    state <- state %*% cpm(monotherapy_ab)
    ab_gr <- expected_gr(state, monotherapy_ab)
    ab_grs <- c(ab_grs, ab_gr)
  }
  AMC_eq_gr <- ab_grs[adaptive_seq_length]
  print(AMC_eq_gr)
  # find adaptive sequence equilibrium growth rate
  adaptive_info <- adaptive_cycling(initial_state, adaptive_seq_length)
  adaptive_grs <- adaptive_info[[2]]
  adaptive_eq_gr <- adaptive_grs[adaptive_seq_length]
  print(adaptive_eq_gr)
  # longterm cost = (AMC monotherapy equilibrium growth rate) - (adaptive sequence equilibrium growth rate)
  longterm_cost <- AMC_eq_gr - adaptive_eq_gr
  
  ## find transient cost
  transient_cost <- 0
  ab_convergence_step <- find_convergence_step(ab_grs)
  adaptive_convergence_step <- find_convergence_step(adaptive_grs)
  convergence_step <- max(ab_convergence_step, adaptive_convergence_step)
  state <- initial_state
  for(n in 1:convergence_step){
    state <- state %*% cpm("AMC")
    AMC_gr <- expected_gr(state, "AMC")
    adaptive_gr <- adaptive_grs[n]
    transient_cost <- transient_cost + (AMC_gr - adaptive_gr)
  }
  return(c(longterm_cost, transient_cost))
}

longterm_costs <- c()
transient_costs <- c()
for(geno in genotypes){
  print(geno)
  initial_state <- rep(0, 16)
  i <- match(geno, genotypes)
  initial_state[i] <- 1
  costs <- find_costs(initial_state, "AMC", adaptive_seq_length = 1000)
  longterm_costs <- c(longterm_costs, costs[1])
  transient_costs <- c(transient_costs, costs[2])
}
names(longterm_costs) <- genotypes
names(transient_costs) <- genotypes

### calculate costs for random initial states###
random_initial_states <- read.table("adaptive_cycling/tables/random_initial_states_1000.txt", sep=",")

rand_is_longterm_costs <- c()
rand_is_transient_costs <- c()
for(i in 1:nrow(random_initial_states)){
  print(i)
  initial_state <- as.matrix(random_initial_states[i,])
  costs_vec <- find_costs(initial_state, "AMC", adaptive_seq_length=1000)
  print("longterm cost")
  print(costs_vec[1])
  print("transient cost")
  print(costs_vec[2])
  rand_is_longterm_costs <-  c(rand_is_longterm_costs, costs_vec[1])
  rand_is_transient_costs <- c(rand_is_transient_costs, costs_vec[2])
}

mean(rand_is_longterm_costs)
sd(rand_is_longterm_costs)

mean(rand_is_transient_costs)
sd(rand_is_transient_costs)