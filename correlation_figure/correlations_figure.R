setwd("~/Antibiotic Cycling")
source(file="cycling.R")
require(graphics)

pdf(file="correlation_figure/ranks_horizontal_logy.pdf",family="Times", width = 7, height = 2,useDingbats=FALSE)

par(mfrow=c(1,3),ps=12)
mar.default <- c(5,4,2,2) + 0.1
par(mar = mar.default + c(0, 2, 0, 0))

single_vs_inf_ranks <- read.table("correlation_figure/single_cycle_vs_inf_cycle_ranks.txt", sep=",")
single_cycle_ranks <- single_vs_inf_ranks[,1]
inf_cycle_ranks <- single_vs_inf_ranks[,2]
#plot(single_cycle_ranks, inf_cycle_ranks, xlab="single cycle\n(short period) rank", ylab=expression(atop("infinite cycle","(short period) rank"%*%10^4)), yaxt="n",ylim=c(-1,300001),pch=16)
plot(single_cycle_ranks, inf_cycle_ranks, xlab="single cycle\n(short period) rank", ylab=expression(atop("infinite cycle","(short period) rank")),pch=16,log="y",ylim=c(10,300000))
#ticks <- c(0,1,5,10,30)
#axis(2, at = ticks * 10000, labels = ticks, las = 1)

lperiod_vs_inf_cycle_ranks <- read.table("correlation_figure/long_period_ranks_vs_inf_cycle_ranks.txt", sep=",")
lperiod_ranks <- lperiod_vs_inf_cycle_ranks[,1]
inf_cycle_ranks <- lperiod_vs_inf_cycle_ranks[,2]
#plot(lperiod_ranks, inf_cycle_ranks,xlab="single cycle\n(long period) rank", ylab=expression(atop("infinite cycle","(short period) rank"%*%10^4)), yaxt="n",ylim=c(-1,300001), pch=16)
plot(lperiod_ranks, inf_cycle_ranks,xlab="single cycle\n(long period) rank", ylab=expression(atop("infinite cycle","(short period) rank")), pch=16,log="y",ylim=c(10,300000))
#ticks <- c(0,10,20,30)
#axis(2, at = ticks * 10000, labels = ticks, las = 1)

speriod_vs_lperiod_ranks <- read.table("correlation_figure/short_period_vs_long_period_ranks.txt", sep=",")
speriod_ranks <- speriod_vs_lperiod_ranks[,1]
lperiod_ranks <- speriod_vs_lperiod_ranks[,2]
#plot(speriod_ranks, lperiod_ranks, xlab="short period\n(single cycle) rank", ylab=expression(atop("long period","(single cycle) rank"%*%10^4)),yaxt="n",ylim=c(-1,300001), pch=16)
plot(speriod_ranks, lperiod_ranks, xlab="short period\n(single cycle) rank", ylab=expression(atop("long period","(single cycle) rank")),log="y",ylim=c(10,300000), pch=16)
#ticks <- c(0,10,20,30)
#axis(2, at = ticks * 10000, labels = ticks, las = 1)

dev.off()

