###
### Example analysis in R:
### Compare the overall survival between the different tumor group
###

library(data.table)
library(survival)

# setwd("Analysis")

# Combine the selected phenotype and survival data into a single table.
data = fread('results/lgg_gbm_survival_input.tsv', data.table=F)
data = data[data$TP53=="WT",]
dim(data)
head(data)
table(data$Expr)

# Create the survival curves for the four different tumor stages.
fit = survfit(Surv(OS_time, OS) ~ Expr, data=data)
r = survdiff(Surv(OS_time, OS) ~ Expr, data=data)
pvalue = pchisq(r$chisq, 3, lower.tail=F)
pvalue

stageColors = c('#e9bb8f', '#8bbad7')
plot(fit, lty=1, col=stageColors, frame.plot=F, lwd=2, xlab='Years', ylab='Survival',main="TCGA LGG/GBM WT p53", mark=-0x00D7L, xaxt='n', yaxt='n', xlim=c(0, 15 * 365))
axis(1, at=seq(0, 15 * 365, by=5 * 365), labels=seq(0, 15, by=5), col='#2f2f2f', col.axis='#2f2f2f')
axis(2, at=seq(0, 1, by=0.2), col='#2f2f2f', col.axis='#2f2f2f', las=1)
text(8 * 365 / 2, 1, paste0('p = ', sprintf('%.4f', pvalue)), adj=c(0.5, 0.5), col='#2f2f2f')
legend('topright', c('High (n=197)', 'Low (n=198)'),col=stageColors, lty=1, lwd=2, bty='n', text.col='#2f2f2f')


