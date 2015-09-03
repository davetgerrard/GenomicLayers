scores <- read.delim("data/HYDRA_runs/auto10k/auto10k.out.tab")
head(scores)

plot(scores); abline(v=3100)
subset(scores, iter > 3000 & iter < 3100)

load("data/HYDRA_runs/auto10k/auto10k.final.Rdata")
length(result)
tail(result$optimScores)

load("data/HYDRA_runs/layer5_10k/currentFactorSet.300.Rdata")
length(currentFactorSet)
#tail(result$optimScores)
