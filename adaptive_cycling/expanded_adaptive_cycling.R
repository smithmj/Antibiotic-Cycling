setwd("~/Antibiotic Cycling") 
source(file="cycling.R")

genotypes = c('0000','1000','0100','0010','0001','1100','1010','1001','0110','0101','0011','1110', '1101','1011','0111','1111')
EX1_growth_rates <- c(.0017, .0017, .001, .0017, .0017, .0017, .0017, .0017, .0017, .0017, .0017, .0017, .001, .0017, .00195, .0017)
EX1_adj_mat <- matrix(0, nrow=16, ncol=16)
rownames(EX1_adj_mat) <- genotypes
colnames(EX1_adj_mat) <- genotypes
# make 0100 a minimum
EX1_adj_mat["0100",] <- c(1,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0)
# make 1101 a minimum
EX1_adj_mat["1101",] <- c(0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1)
# make 0111 a maximum
EX1_adj_mat["1111","0111"] <- 1
EX1_adj_mat["0011","0111"] <- 1
EX1_adj_mat["0101","0111"] <- 1
EX1_adj_mat["0110","0111"] <- 1
# if a row is all 0, replace the diagonal entery with a 1
for(i in 1:nrow(EX1_adj_mat)){
  if(sum(EX1_adj_mat[i,]) == 0){
    EX1_adj_mat[i,i] <- 1
  }
}

expanded_antibiotics <- c(antibiotics, "EX1")
expanded_growth_rates <- rbind(growth_rates, EX1_growth_rates)
rownames(expanded_growth_rates) = expanded_antibiotics

initial_state <- c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#adaptive_optimization(initial_state,ab_list=expanded_antibiotics,gr=expanded_growth_rates) 
adaptive_list <- adaptive_cycling(initial_state, 300, ab_list=expanded_antibiotics, gr=expanded_growth_rates)
adaptive_seq <- adaptive_list[[1]]