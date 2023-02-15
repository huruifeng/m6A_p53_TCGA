###
### Example analysis in R:
### Compare the overall survival between the different tumor group
###



library(data.table)
library(survival)


# Combine the selected phenotype and survival data into a single table.
data = fread('survival_input21.tsv', data.table=F)
dim(data)
head(data)
table(data$Group)

# Create the survival curves for the four different tumor stages.
fit = survfit(Surv(OS.time, OS) ~ Group, data=data)
r = survdiff(Surv(OS.time, OS) ~ Group, data=data)
pvalue = pchisq(r$chisq, 3, lower.tail=F)
pvalue

stageColors = c('#FF0000', '#0000FF')
plot(fit, lty=1, col=stageColors, frame.plot=F, lwd=3, xlab='Days', ylab='Survival',mark.time=TRUE, pch=3, xaxt='n', yaxt='n')
axis(1,  col='#2f2f2f', col.axis='#2f2f2f')
axis(2, at=seq(0, 1, by=0.2), col='#2f2f2f', col.axis='#2f2f2f', las=1)
# text(8 * 365 / 2, 1, paste0('p = ', sprintf('%3.2e', pvalue)), adj=c(0.5, 0.5), col='#2f2f2f')
title("Kaplan-Meier Curve - LGG/GBM")
# legend('topright', c('High (n=345)', 'Low (n=346)'),col=stageColors, lty=1, lwd=2, bty='n', text.col='#2f2f2f')


