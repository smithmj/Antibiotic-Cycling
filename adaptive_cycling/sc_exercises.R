#####################################################################
# Exercises to investigate potential advantages of cycling
#####################################################################

#setwd("~/Antibiotic Cycling")
require(stringdist)
source("cycling.R")

# Create adjacency matrix given growth rates
# Note: saved to file for use with cpm() function
make_adj <- function(ab_name,ab_grs,gts=genotypes){
	A = matrix(0,nrow=length(gts),ncol=length(gts))
	rownames(A) = gts
	colnames(A) = gts
	for (i in gts) {
		for (j in gts) {
			if (stringdist(i,j,method="hamming")==1 && ab_grs[i] < ab_grs[j]) {
		        	A[i,j]=1 }
		    }
		    if (sum(A[i,])==0) {	
		    	A[i,i]=1 }
		    }
		    filename <- paste(ab_name,'_adj_mat.csv',sep="")
		    write.table(A, file = filename,row.names=FALSE,col.names=FALSE, sep=",")
# filename <- paste(ab_name,'_growth_rates.csv',sep="")
# write.table(t(fab_growth_rates),file=filename,row.names=FALSE,col.names=FALSE,sep=",")
return(A)
}

# Add new drug info to antibiotics list and growth rates
# (changes are made to global workspace variables)
add_drug <- function(ab_name,ab_grs) {
	antibiotics <<- c(antibiotics,ab_name)
	gr_names <- row.names(growth_rates)
	gr_names <- c(gr_names,ab_name)
	growth_rates <<- rbind(growth_rates,ab_grs)
	rownames(growth_rates) <<- gr_names
}

# Save for easy resetting
original_antibiotics <- antibiotics

#####################################################################
# Can we obtain cycling with two drugs (with distant maxima)?
# Answer: Yes

# Part 1: Create new antibiotic landscape from existing one
# by reversing AMC landscape

drugname <- 'X1'
index <- match('AMC',antibiotics) 
fab_growth_rates <- growth_rates[index,]
sorted_gr <- sort(fab_growth_rates)
for (a in 1:length(sorted_gr)) {
	g1 = names(sorted_gr)[a]
	g2 = names(sorted_gr)[length(sorted_gr)-a+1]
	fab_growth_rates[g1] <- sorted_gr[g2]
}

make_adj(drugname,fab_growth_rates) 
# Move new file to adj_mat dir
add_drug(drugname,fab_growth_rates)
inits <- runif(length(genotypes)); inits <- inits/sum(inits)
selected_cycles <- adaptive_cycling(inits,100)

# > selected_cycles
# [[1]]
# [1] "AMC" "X1"  "AMC" "X1"  "AMC" "X1"  "AMC" "X1"  "AMC" "X1"

# [[2]]
# [1] "0.00166426448336328" "0.00160520405748432" "0.00165589702174815" "0.00160893248390516"
# [5] "0.00165742164483493" "0.00160990814422393" "0.00165812798402371" "0.00161076629585316"
# [9] "0.00165824649851199" "0.00161145778782453"

# [[3]]
# 0000       1000 0100       0010        0001      1100       1010      1001     0110
# [1,] 0.04504897 0.02342961    0 0.07683229 0.005852438 0.1394732 0.01121221 0.1865638 0.117407
# 0101       0011       1110 1101       1011      0111         1111
# [1,] 0.09559804 0.02543395 0.08176768    0 0.03801648 0.1530036 0.0003607498


# Part 2: New antibiotic is like previous but has slightly worse
# maximum growth rate instead
growth_rates <- growth_rates[original_antibiotics,]
drugname <- 'X2'

# ...as before
index <- match('AMC',antibiotics) 
fab_growth_rates <- growth_rates[index,]
sorted_gr <- sort(fab_growth_rates)
for (a in 1:length(sorted_gr)) {
	g1 = names(sorted_gr)[a]
	g2 = names(sorted_gr)[length(sorted_gr)-a+1]
	fab_growth_rates[g1] <- sorted_gr[g2]
}

# > max(growth_rates['AMC',])
# [1] 0.001914167

idx = which(fab_growth_rates==max(growth_rates['AMC',]))
fab_growth_rates[idx] <- fab_growth_rates[idx]*1.1
make_adj(drugname,fab_growth_rates) 
add_drug(drugname,fab_growth_rates)
inits <- runif(length(genotypes)); inits <- inits/sum(inits)
selected_cycles <- adaptive_cycling(inits,10)

# > selected_cycles
# [[1]]
 # [1] "AMC" "X2"  "CTX" "AMC" "X2"  "AMC" "X2"  "AMC" "X2"  "AMC"

# [[2]]
 # [1] "0.00168172432418053" "0.00159304616335581" "0.0015670569440571"  "0.00161973234013948"
 # [5] "0.00162806355944335" "0.00162751526036906" "0.00164130644644135" "0.00164129116919903"
 # [9] "0.00164452544519165" "0.00164896645492135"

# [[3]]
           # 0000       1000      0100       0010       0001        1100       1010 1001
# [1,] 0.02297303 0.01707503 0.1510865 0.03849572 0.06799768 0.007613535 0.06983776    0
             # 0110       0101       0011       1110      1101       1011 0111      1111
# [1,] 0.0001888511 0.05472521 0.06645839 0.04807875 0.2533778 0.06637856    0 0.1357131


# What is effect on equilibrium transition matrix?
X2_tm <- cpm('X2')
AMC_tm <- cpm('AMC')
cycle_tm <- AMC_tm %*% X2_tm
cycle_eq_tm <- cycle_tm %^% 100000

heatmap(cycle_eq_tm,Rowv=NA,Colv=NA,revC=TRUE)

# > max(cycle_eq_tm)
# [1] 0.458842

AMC_eq_tm <- AMC_tm %^% 100000
heatmap(AMC_eq_tm,Rowv=NA,Colv=NA,revC=TRUE)
AMC_eq_states <- inits%*%AMC_eq_tm
expected_gr(AMC_eq_states,'AMC')

# > expected_gr(AMC_eq_states,'AMC')
            # [,1]
# [1,] 0.001863666

cycle_eq_states <- inits%*%cycle_eq_tm  # averaging states?
expected_gr(cycle_eq_states,'AMC')

# > expected_gr(cycle_eq_states,'AMC')
#            [,1]
# [1,] 0.001041613

# > expected_gr(cycle_eq_states,'X2')
            # [,1]
# [1,] 0.001624786


