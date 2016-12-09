#setwd("~/Antibiotic Cycling")
library(expm)
library(viridis)

cutoff <- 0.001 # Growth rate at or above which better antibiotic sought
antibiotics = c('AMP','AM','CEC','CTX','ZOX','CXM','CRO','AMC', 'CAZ','CTT','SAM','CPR','CPD','TZP','FEP')
genotypes = c('0000','1000','0100','0010','0001','1100','1010','1001','0110','0101','0011','1110', '1101','1011','0111','1111')

growth_rates = as.matrix(read.csv("growth_rate.csv", sep = ",", header = FALSE))
rownames(growth_rates) = antibiotics
colnames(growth_rates) = genotypes

## Imports the adjacency matrix of a given antibiotic from the adj_mat directory ##
import_adj_mat <- function(ab){
  file_name = paste("adj_mat/", ab, "_adj_mat.csv", sep = "")
  mat <- as.matrix(read.csv(file_name, sep=",", header = F))
  rownames(mat) = genotypes
  colnames(mat) = genotypes
  return(mat)
}

## Creates the matrix of transition probabilities for a given adjacency matrix ##
## according to the correlated probability model                               ##
## --- as of now, this does not have immigration probability --- ##
cpm <- function(ab, gr=growth_rates){
	#  print(gr)
  adj_mat <- import_adj_mat(ab)
  subt_mat <- matrix(0, nrow=ncol(adj_mat), ncol=ncol(adj_mat))
  for(u in 1:nrow(subt_mat)){
    for(v in 1:ncol(subt_mat)){
		#      print(gr[ab,v])
      if(adj_mat[u,v] == 1 && gr[ab, v] > gr[ab, u]){
        num = gr[ab, v] - gr[ab, u]
        denom = 0
        for(j in 1:ncol(adj_mat)){
          if(adj_mat[u,j] == 1 && gr[ab,j] > gr[ab,u]){
            denom = denom + gr[ab,j] - gr[ab,u]
          }
        }
        subt_mat[u,v] = num / denom
      }
    }
    if(sum(subt_mat[u,]) == 0){
      subt_mat[u,u] = 1
    }
  }
  rownames(subt_mat) <- genotypes
  colnames(subt_mat) <- genotypes
  return(subt_mat)
}

cycle_mat <- function(ab_seq, long_period=F){
  M <- diag(length(genotypes))
  for(ab in ab_seq){
     if(long_period){
      M_ab <- eq_mat(cpm(ab))
    } else {
      M_ab <- cpm(ab)
    }
    M <- M %*% M_ab
  }
  return(M)
}

## Finds the equilibrium transition matrix given a transition probability matrix ##
eq_mat <- function(mat){
  n <- 100000
  return(mat %^% n)
}

## Finds the expected growth rate given the genotype state and the antibiotic ##
expected_gr <- function(state, ab, gr=growth_rates){
  ex_gr <- state %*% gr[ab,]
  return(ex_gr)
}

adaptive_optimization <- function(state,gr_cutoff=cutoff,ab_list=antibiotics,gr=growth_rates){
  min_gr <- Inf
  best_ab <- NA
  gr_vec <- c()
  for(ab in ab_list){
    new_state <- state %*% cpm(ab,gr=gr)
    ab_gr <- expected_gr(new_state,ab,gr=gr)
    gr_vec <- c(gr_vec, ab_gr)
    if(ab_gr <= min_gr){
      best_ab <- ab
      min_gr <- ab_gr 
    }  
  }
  ## print a warning if there was a tie ##
  sorted_grs <- sort(gr_vec)
  if(sorted_grs[1] == sorted_grs[2]){
    print("The minimum growth rate is the same for two or more drugs")
  }
  
  return(c(best_ab, min_gr))
}

adaptive_cycling <- function(initial_state, num_cycles, gr_cutoff=cutoff, ab_list=antibiotics,gr=growth_rates){
  adaptive_seq <- c()
  seq_growth_rates <- c()
  while(length(adaptive_seq) < num_cycles){
    opt_info <- adaptive_optimization(initial_state,gr_cutoff=gr_cutoff,ab_list=ab_list,gr=gr)
    best_ab <- opt_info[1]
    adaptive_seq <- c(adaptive_seq, best_ab)
    ex_gr <- opt_info[2]
    seq_growth_rates <- c(seq_growth_rates, ex_gr)
    initial_state <- initial_state %*% cpm(best_ab, gr=gr)
    while(ex_gr <= gr_cutoff && length(adaptive_seq) < num_cycles){
      adaptive_seq <- c(adaptive_seq, best_ab)
      initial_state <- initial_state %*% cpm(best_ab)
      ex_gr <- expected_gr(initial_state, best_ab, gr=gr)
      seq_growth_rates <- c(seq_growth_rates, ex_gr)
    }
  }
  adaptive_list <- list(adaptive_seq, as.numeric(seq_growth_rates), initial_state)
  return(adaptive_list)
}

plot_cycles <- function(ab_seq,initial_state,combine_cycle=F,long_period=F,num_cycles = length(ab_seq),xlab="",ylab="frequency",xticklab=NA,legend=T){
  
  if(combine_cycle){
    #in this case, we want each bar on the graph to represent the outcome from 
    #from the entire cycle
    plot_info <- matrix(0,nrow=length(genotypes),ncol=num_cycles+1)
    plot_info[,1] <- initial_state
    mat_cycle <- cycle_mat(ab_seq,long_period)
    for(i in 1:num_cycles){
      initial_state <- initial_state %*% mat_cycle
      plot_info[,i+1] <- initial_state      
    }
  } else {
    #in this case we want each bar on the graph to represent the outcome from
    #a single antibiotic treatment
    plot_info <- matrix(0,nrow=length(genotypes),ncol=length(ab_seq)+1)
    plot_info[,1] <- initial_state
    for(i in 1:length(ab_seq)){
      ab <- ab_seq[i]
      if(long_period){
        M_ab <- eq_cycle(cpm(ab))
      } else {
        M_ab <- cpm(ab)
      }
      initial_state <- initial_state %*% M_ab
      plot_info[,i + 1] <- initial_state
    }
  }
  
  #colors <- viridis(length(genotypes))
  colors <- c("#B550A8", "#BE4F8E", "#BD556E", "#B45F47", "#A66A00", "#917400", "#767E00", "#508500", "#008B34", "#008F5E", "#00907F", "#008D9B", "#0086B1", "#2F7BBF", "#786BC2", "#9E5BBA")
  barplot(plot_info,col=colors,space=c(.1,.1),xlab=xlab,ylab=ylab,las=1)
  if(legend){
    legend("topright", inset=c(-0.06,-.05), genotypes, bty="n", fill=colors,xpd=NA)
  }
  if(is.na(xticklab)){
    ticks <- seq(0,num_cycles,by=2)
    tick_loc <- ticks + (ticks * 0.1 + 0.1) + 0.5
  } else {
    ticks <- xticklab
    tick_loc <- seq(0,num_cycles)
    tick_loc <- tick_loc + (tick_loc * 0.1 + 0.1) + 0.5
  }
  axis(side=1,labels=ticks,at=tick_loc)
}

plot_growth_rates <- function(ab_seq, initial_state){
  growth_rates <- c()
  for(ab in ab_seq){
    gr <- expected_gr(initial_state, ab)
    growth_rates <- c(growth_rates, gr)
    initial_state <- initial_state %*% cpm(ab)
  }
  plot(growth_rates, pch=ncol(growth_rates), type="o", xlim=c(0,length(ab_seq)), xlab="", xaxt="n", yaxt="n", ylab=expression("expected growth rate"%*%10^-3), ylim=c(0,0.002))
  ticks <- c(0.0, 0.5, 1.0, 1.5, 2.0)
  axis(2, at = ticks / 1000, labels = format(ticks, format="f", digits=2))
}

find_max_cost <- function(ab, initial_geno, onset_step){
  initial_state <- rep(0, length(genotypes))
  #find the index of the initial_geno in the vector genotypes.
  #this will give the position of the initial state that should
  #be a 1
  index <- match(initial_geno, genotypes)
  initial_state[index] <- 1
  adaptive_list <- adaptive_cycling(initial_state, onset_step)
  adaptive_grs <- as.numeric(adaptive_list[[2]])
  print("adaptive_grs")
  print(adaptive_grs)
  fixed_grs <- c()
  for(i in 1:onset_step){
    initial_state <- initial_state %*% cpm(ab)
    gr <- expected_gr(initial_state, ab)
    fixed_grs <- c(fixed_grs, gr)
  }
  diff_grs <- fixed_grs - adaptive_grs
  return(max(diff_grs))
}
