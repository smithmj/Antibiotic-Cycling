setwd("~/Antibiotic Cycling")
source(file="cycling.R")


pdf(file="short_period_vs_long_period/short_period_vs_long_period.pdf",family="Times", width = 7.28, height = 4, pointsize=9)

par(mfrow=c(1,2))
initial_state = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
par(mai=c(1.02,.8,0.7,0))
plot_cycles(c("AM", "TZP"),initial_state,xlab="AM-TZP short period cycle",num_cycles=15,combine_cycle=T,legend=F)
par(mai=c(1.02,.3,0.7,.5))
plot_cycles(c("AM","TZP"),initial_state,long_period=T,xlab="AM-TZP long period cycle",num_cycles=15,combine_cycle=T)

dev.off()
