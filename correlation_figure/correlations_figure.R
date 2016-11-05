setwd("~/Antibiotic Cycling")
source(file="cycling.R")
require(graphics)

pdf(file="correlation_figure/ranks_vertical_labeled.pdf",family="Times", width = 2.5, height = 7)

par(mfrow=c(3,1))
mar.default <- c(5,4,2,2) + 0.1
par(mar = mar.default + c(0, 2, 0, 0)) 

single_vs_inf_ranks <- read.table("correlation_figure/single_cycle_vs_inf_cycle_ranks.txt", sep=",")
single_cycle_ranks <- single_vs_inf_ranks[,1]
inf_cycle_ranks <- single_vs_inf_ranks[,2]
plot(single_cycle_ranks, inf_cycle_ranks, xlab="single cycle\n(short period) rank", ylab=expression(atop("infinite cycle","(short period) rank"%*%10^4)), yaxt="n",pch=16)
ticks <- c(0,5,10,15,20,25,30)
axis(2, at = ticks * 10000, labels = ticks, las = 1)

lperiod_vs_inf_cycle_ranks <- read.table("correlation_figure/long_period_ranks_vs_inf_cycle_ranks.txt", sep=",")
lperiod_ranks <- lperiod_vs_inf_cycle_ranks[,1]
inf_cycle_ranks <- lperiod_vs_inf_cycle_ranks[,2]
plot(lperiod_ranks, inf_cycle_ranks,xlab="single cycle\n(long period) rank", ylab=expression(atop("infinite cycle","(short period) rank"%*%10^4)), yaxt="n",pch=16)
ticks <- c(0,5,10,15,20,25,30)
axis(2, at = ticks * 10000, labels = ticks, las = 1)

speriod_vs_lperiod_ranks <- read.table("correlation_figure/short_period_vs_long_period_ranks.txt", sep=",")
speriod_ranks <- speriod_vs_lperiod_ranks[,1]
lperiod_ranks <- speriod_vs_lperiod_ranks[,2]
plot(speriod_ranks, lperiod_ranks, xlab="short period\n(single cycle) rank", ylab=expression(atop("long period","(single cycle) rank"%*%10^4)),yaxt="n",pch=16)
ticks <- c(0,5,10,15,20,25,30)
axis(2, at = ticks * 10000, labels = ticks, las = 1)

dev.off()

