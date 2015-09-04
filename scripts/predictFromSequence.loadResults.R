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


load("data//HYDRA_runs//layer5_10k//layer5_10k.final.Rdata")
length(result)
table.factorSet(factorSetRandom)
type.patterns(table.factorSet(factorSetRandom))
table.factorSet(result[1:30])
type.patterns(table.factorSet(result[1:30]))

head(result$optimScores)
plot(result$optimScores); abline(v=160); abline(v=170)
table.factorSet(result[1:30])  # looking for factors added between 160 and 170, although tab file suggests between 200-300
