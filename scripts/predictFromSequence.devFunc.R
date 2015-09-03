## Scratchbox to test the latest version of the functions.


require(Biostrings)
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')
source('scripts/predictFromSequence.functions.R')





## Try some alternative algorithms for running a series of mods over a sequence.
# Main target is improved speed.

targetSeq <- readDNAStringSet("C:/Users/Dave/Dropbox/Temp/hg19.HOXA.fa")[[1]]
emptyLayer <-  BString(paste(rep(0,length(targetSeq)), collapse=""))



n.layers <- 1

layerSet.1 <- list(LAYER.0 = targetSeq)
for(i in 1:n.layers) {
  layerSet.1[[paste('LAYER.', i , sep="")]] <- emptyLayer
}

layerList.1 <- list(layerSet=layerSet.1, history=NULL)

n.factors <- 30

#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), n.factors, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T)
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T, prob=c(10,10,2,2))
factorSetRandom <- list()
for(i in 1:n.factors) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), layerSet.1, type=bindingFactorTypes[i], test.layer0.binding=FALSE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)


system.time(modLayerSet <- runLayerBinding(layerList=layerList.1, factorSet = factorSetRandom, iterations=100))   # 3 mins.
lapply(modLayerSet$layerSet[-1], letterFrequency, letters="1")    # 3013 mods

system.time(modLayerSet.fast <- runLayerBinding.fast(layerList=layerList.1, factorSet = factorSetRandom, iterations=100))   # 1 min
lapply(modLayerSet.fast$layerSet[-1], letterFrequency, letters="1")  # 4621 mods

system.time(modLayerSet.fast10k <- runLayerBinding.fast(layerList=layerList.1, factorSet = factorSetRandom, iterations=10000))   # 30 secs!  
lapply(modLayerSet.fast10k$layerSet[-1], letterFrequency, letters="1")  # 125,239 mods

system.time(modLayerSet.fast10k.2 <- runLayerBinding.fast(layerList=layerList.1, factorSet = factorSetRandom, iterations=10000))   # 30 secs!  
lapply(modLayerSet.fast10k.2$layerSet[-1], letterFrequency, letters="1")  # 123, 450 mods

# check correlation between the first layers of each of  the reps.
cor.layers( modLayerSet.fast10k.2$layerSet[['LAYER.1']],  modLayerSet.fast10k$layerSet[['LAYER.1']])




# see if optimiseFactorSet will work with new fast function

transcriptTable <- read.delim("data/hg19.HOXA.transcripts.tab")
head(transcriptTable)
base.0 <- 27106375   # one before the start co-ordinate used to get this DNA sequence.

tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.

tss.vector <- rep(0, nchar(targetSeq))
tss.vector[tss.positions] <- 1

# widen the tss positions to simulate promoters?
upstream.prom <- 200
downstream.prom <- 200
for(i in 1:nrow(transcriptTable)) {
  tss.position <- ifelse(transcriptTable$strand[i] == "+", transcriptTable$txStart[i], transcriptTable$txEnd[i]) - base.0
  prom.start <- ifelse(transcriptTable$strand[i] == "+", tss.position - upstream.prom, tss.position - downstream.prom)
  prom.end <- ifelse(transcriptTable$strand[i] == "-", tss.position + upstream.prom, tss.position + downstream.prom)
  tss.vector[prom.start:prom.end] <- 1
}
sum(tss.vector)
length(tss.vector)


test_function <- function(layerList, targetLayer="Layer.1", target.vec)  {
  layer.vec <- as.numeric(strsplit(as.character(layerList$layerSet[[targetLayer]]),"")[[1]])
  return(cor(target.vec, layer.vec))
}


n.iter <- 1000
system.time(result <- optimiseFactorSet(layerList=layerList.1, factorSetRandom, testing.function=test_function, target.layer="LAYER.1", target.vec=tss.vector, n.iter=n.iter, mut.rate=0.1, modsPerCycle=10000))

plot(result$optimScores)

save(factorSetRandom, result, file=paste("results/run.fast", n.iter,"Rdata", sep="."))


modLayerList <- runLayerBinding.fast(layerList=layerList.1, factorSet = result[-(length(result))], iterations=10000)  # result[-(length(result))] because last object is a set of scores
test_function(modLayerList,  targetLayer="LAYER.1", tss.vector)  # still quite variable


nrow(transcriptTable)  # 72 TSS in test set. (though some are overlapping)
length(tss.positions)  # 71 unique

targetLayer <- "LAYER.1"
result.vector <- as.numeric(strsplit(as.character(modLayerList$layerSet[[targetLayer]]),"")[[1]])
plot(tss.vector)
points(result.vector-.05, col="blue")

mods.list <- list()
for(i in 1:10) {
  mods.list[[i]] <- runLayerBinding.fast(layerList=layerList.1, factorSet = result[-(length(result))], iterations=10000)  # result[-(length(result))] because last object is a set of scores
  
}

plot(tss.vector, xlim=c(50000, 60000))
for(i in 1:10) {
  result.vector <- as.numeric(strsplit(as.character(mods.list[[i]]$layerSet[[targetLayer]]),"")[[1]])
  points(result.vector-(.05*i), col="blue")
}

# or the whole thing
plot(tss.vector)
for(i in 1:10) {
  result.vector <- as.numeric(strsplit(as.character(mods.list[[i]]$layerSet[[targetLayer]]),"")[[1]])
  points(result.vector-(.05*i), col="blue")
}

# TODO
# take the mean value for each bp across the mods.list
# plot in a sub-region

# use more, don't keep the full result calc sum along the way.
reps <- 100
result.vector <- rep(0, length(tss.vector))
for(i in 1:reps) {
  print(".")
  thisModLayerList <- runLayerBinding.fast(layerList=layerList.1, factorSet = result[-(length(result))], iterations=10000)  # result[-(length(result))] because last object is a set of scores
  this.result.vector <- as.numeric(strsplit(as.character(thisModLayerList$layerSet[[targetLayer]]),"")[[1]])
  result.vector <- result.vector + this.result.vector
}

plot(tss.vector* 1.05, xlim=c(50000, 60000))
points(result.vector/reps, col="blue")
points(which(tss.vector == 1), (result.vector/reps)[tss.vector == 1], col="green")


plot(tss.vector* 1.05)
points(result.vector/reps, col="blue")
points(which(tss.vector == 1), (result.vector/reps)[tss.vector == 1], col="green")

plot(tss.vector* 1.05, xlim=c(100000, 120000))
points(result.vector/reps, col="blue")
points(which(tss.vector == 1), (result.vector/reps)[tss.vector == 1], col="green")

plot(which(tss.vector == 1), rep(1.05, sum(tss.vector)), xlim=c(80000, 150000), ylim=c(0,1.1))
points(result.vector/reps, col="blue")
points(which(tss.vector == 1), (result.vector/reps)[tss.vector == 1], col="green")

png(paste("figures/run.fast", n.iter,"png", sep="."), width=4000, height=1000, res=150)
plot(which(tss.vector == 1), rep(1.05, sum(tss.vector)), xlim=c(0, length(tss.vector)), ylim=c(0,1.1), xlab="position", ylab="Proportion of mods base is altered")
points(result.vector/reps, col="blue")
points(which(tss.vector == 1), (result.vector/reps)[tss.vector == 1], col="green")
dev.off()

# can use the summed values to re-score the factors
sum(result.vector/reps > .4)
sum(result.vector/reps)
sum(tss.vector)
cor(tss.vector, result.vector/reps)

# does a threshold improve the correlation?  - not much.
threshold <- .4
top.vector <- rep(0, length(tss.vector))
top.vector[(result.vector/reps) > threshold] <- 1
cor(tss.vector, top.vector)
cor.scores <- numeric()
threshes <- seq(0,1,0.1)
for(i in 1:length(threshes)) {
  top.vector <- rep(0, length(tss.vector))
  top.vector[(result.vector/reps) > threshes[i]] <- 1
  cor.scores[i] <- cor(tss.vector, top.vector)
}

pdf(paste("figures/run.fast", n.iter,"correlation.threshold.pdf", sep="."), width=8, height=8)
plot(threshes,cor.scores, xlab="Threshold for proportion of runs picked", ylab="Correlation")
abline(h=cor(tss.vector, result.vector/reps), lty=2)
dev.off()

print.bfSet(factorSetRandom)

# test that logfiles work

n.iter <- 20
system.time(testResult <- optimiseFactorSet(layerList=layerList.1, factorSetRandom, testing.function=test_function, target.layer="LAYER.1", target.vec=tss.vector, n.iter=n.iter, mut.rate=0.1, modsPerCycle=10000, logFile="results/test.log", logCycle=5))

# test whether 'maxNoChange' working
system.time(testResult <- optimiseFactorSet(layerList=layerList.1, factorSetRandom, testing.function=test_function, target.layer="LAYER.1", target.vec=tss.vector, n.iter=n.iter, mut.rate=0.1, modsPerCycle=10000, logFile="results/test.log", logCycle=5,maxNoChange=2))


########## need better characterisation of the resulting factor set
# Would like tale of which layers are targets and which are mods.
# Many to may relationship - each factor could recognise multiple layers and modify multiple.
# Maybe a eulergrid, scoring each factor for both its motif targets and mod targets (and it's type?).


factorSet <- factorSetRandom



type.table <- table.factorSet(factorSetRandom)


table(type.table$group, type.table$layer)



type.patterns(type.table)
type.patterns(table.factorSet(factorSetRandom))
type.patterns(table.factorSet(result[1:30]))
type.patterns(table.factorSet(currentFactorSet))


# plot a euler grid of the mod and target combinations to see what they look like.
source("c:/users/dave/utilsGerrardDT/dataToEulerGrid.R")
type.pattern.table <- as.data.frame(type.patterns(table.factorSet(currentFactorSet)))
type.pattern.table <- type.pattern.table[, sort(names(type.pattern.table))]
euler.table <- scoreCardinalities(type.pattern.table)
par(mar=c(3,10,3,3))
plotEuler(subset(euler.table, select=- count), counts=euler.table$count)

