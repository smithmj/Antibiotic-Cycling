setwd("~/Antibiotic Cycling")
source(file="cycling.R")

pdf(file="strategy_summary_figure/strategy_summary_figure.pdf",family="Times", width = 3.54, height = 3.54,pointsize=9)

## need to write the code that geneartes these numbers ##
no_imm <- c(.70, .60, .15, .14)
imm <- c(.70, .57, .17, .11)
mat <- matrix(c(no_imm, imm), nrow = 2, byrow = T)

barpos <- barplot(mat, beside = T, xlim=c(.5,9.5), space=c(0,0.5), col=c("gray22", "gray62"),ylab="final frequency sensitive", ylim=c(0,.8),names.arg=c("s-c\ns-p","i-c\ns-p","s-c\nl-p","i-c\nl-p"), las=1)
text(barpos,mat+.05,labels=formatC(mat,format="f", digits=2)) 
legend("topright", c("no immigration","immigration"), bty="n", fill=c("gray22", "gray62"))

dev.off()
