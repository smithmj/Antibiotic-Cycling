setwd("~/Antibiotic Cycling")
source(file="cycling.R")


pdf(file="cycle_examples/cycle_examples.pdf",family="Times", width = 7, height = 3, pointsize=9)

par(mfrow=c(1,2), las = 1)
initial_state = c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#bottom,left,top,right
par(mai=c(1,0.6,0.1,0.7))
plot_cycles(c("AM", "TZP"),initial_state,xlab="AM-TZP cycles",num_cycles=15,combine_cycle=T,legend=F)
par(mai=c(1,0.6,0.1,0.7))
plot_cycles(c("CEC","SAM","CTX","AM"),initial_state,xlab="CEC-SAM-CTX-AM cycles",num_cycles=15,combine_cycle=T,legend=F)

colors <- c("#B550A8", "#BE4F8E", "#BD556E", "#B45F47", "#A66A00", "#917400", "#767E00", "#508500", "#008B34", "#008F5E", "#00907F", "#008D9B", "#0086B1", "#2F7BBF", "#786BC2", "#9E5BBA")
legend("right", inset=c(-0.3,0),genotypes, bty="n", fill=colors, horiz=F, xpd=NA)


dev.off()
