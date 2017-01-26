setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')


model.methods <- c("NNN", "CpG", "GpC")

scores.null <- scan("results/mus_X_inactivation/mXi.2016.12.09.null.CpG/scores.30.txt")

scores.CpG <- scan("results/mus_X_inactivation/mXi.2016.12.09.basic.CpG/scores.30.txt")

scores.GpC <- scan("results/mus_X_inactivation/mXi.2016.12.09.invert.GpC/scores.30.txt")


png(file="results/mus_X_inactivation/mXi.2016-12-09.compareModels.bxp.png", width=600, height=800, res=150)
boxplot(scores.null, xlim=c(0,4), ylim=c(-1,1)) ; mtext("NNN", side=1, at =1, line=2)
boxplot(scores.CpG, add=T, at=2); mtext("CpG", side=1, at =2, line=2)
boxplot(scores.GpC, add=T, at=3); mtext("GpC", side=1, at =3, line=2)
abline(h=c(-1, -.5, 0, .5, 1), lty=2, col="grey")
dev.off()

wilcox.test(scores.CpG, scores.null, alternative = "greater")
wilcox.test(scores.CpG, scores.GpC, alternative = "greater")
