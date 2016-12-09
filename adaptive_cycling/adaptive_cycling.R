setwd("~/Antibiotic Cycling")
source(file="cycling.R")

pdf(file="adaptive_cycling/adaptive_cycling_examples.pdf",family="Times", width = 7.28, height = 4.5, pointsize=9,useDingbats=FALSE)

par(mfrow=c(2,3), las=2)
par(las=2)

initial_state_1 <- c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
adaptive_list <- adaptive_cycling(initial_state_1, 20)
adaptive_seq_1 <- adaptive_list[[1]]
ticks_1 <- c("", adaptive_seq_1)
par(mai=c(0.1,0.5,0.82,0.2))
plot_growth_rates(adaptive_seq_1, initial_state_1)

initial_state_2 <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0)
adaptive_list <- adaptive_cycling(initial_state_2, 20)
adaptive_seq_2 <- adaptive_list[[1]]
ticks_2 <- c("", adaptive_seq_2)
par(mai=c(0.1,0.5,0.82,0.2))
plot_growth_rates(adaptive_seq_2, initial_state_2)

initial_state_3 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
adaptive_list <- adaptive_cycling(initial_state_3, 20)
adaptive_seq_3 <- adaptive_list[[1]]
ticks_3 <- c("", adaptive_seq_3)
par(mai=c(0.1,0.5,0.82,0.2))
plot_growth_rates(adaptive_seq_3, initial_state_3)

colors <- c("#B550A8", "#BE4F8E", "#BD556E", "#B45F47", "#A66A00", "#917400", "#767E00", "#508500", "#008B34", "#008F5E", "#00907F", "#008D9B", "#0086B1", "#2F7BBF", "#786BC2", "#9E5BBA")
par(mai=c(1.02,0.5,0.1,0.2))
plot_cycles(adaptive_seq_1, initial_state_1, xlab="step", xticklab=ticks_1, legend=F)
par(mai=c(1.02,0.5,0.1,0.2))
plot_cycles(adaptive_seq_2, initial_state_2, xlab="step", xticklab=ticks_2, legend=F)
legend("bottom", inset=c(0,-.6), genotypes, bty="n", fill=colors, horiz=T, xpd=NA)
par(mai=c(1.02,0.5,0.1,0.2))
plot_cycles(adaptive_seq_3, initial_state_3, xlab="step", xticklab=ticks_3, legend=F)
#colors <- jet.colors(16)
#plot.new()
par(xpd=NA)

dev.off()


