setwd("~/Antibiotic Cycling")
source(file="cycling.R")


pdf(file="cycle_examples/cycle_examples.pdf",family="Times", width = 7.28, height = 4, pointsize=9)

par(mfrow=c(1,2), las = 1)
initial_state = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#bottom,left,top,right
par(mai=c(1.02,.8,0.7,0))
plot_cycles(c("AM", "TZP"),initial_state,xlab="AM-TZP cycles",num_cycles=15,combine_cycle=T,legend=F)
par(mai=c(1.02,.3,0.7,.5))
plot_cycles(c("CEC","SAM","CTX","AM"),initial_state,xlab="CEC-SAM-CTX-AM cycles",num_cycles=15,combine_cycle=T)

dev.off()
