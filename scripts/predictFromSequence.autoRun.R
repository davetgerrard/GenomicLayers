

# Keep running optimisation until some threshold is passed or until a set number of runs have passed without improvement
# Log results every x iterations.



require(Biostrings)
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')
source('scripts/predictFromSequence.functions.R')





## Try some alternative algorithms for running a series of mods over a sequence.
# Main target is improved speed.

targetSeq <- readDNAStringSet("C:/Users/Dave/Dropbox/Temp/hg19.HOXA.fa")[[1]]
transcriptTable <- read.delim("data/hg19.HOXA.transcripts.tab")
base.0 <- 27106375 
n.layers <- 1
n.factors <- 30
upstream.prom <- 200
downstream.prom <- 200
n.iter <- 1000
mut.rate <- 0.1
modsPerCycle <- 10000
logCycle<- 100 
maxNoChange<- 1000
outPutDir <- "."
runName <- "test"
logFile <- paste(outputDir, "/" , runName, ".out.tab", sep="")



emptyLayer <-  BString(paste(rep(0,length(targetSeq)), collapse=""))



layerSet.1 <- list(LAYER.0 = targetSeq)
for(i in 1:n.layers) {
  layerSet.1[[paste('LAYER.', i , sep="")]] <- emptyLayer
}

layerList.1 <- list(layerSet=layerSet.1, history=NULL)



#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), n.factors, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T)
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T, prob=c(10,10,2,2))
factorSetRandom <- list()
for(i in 1:n.factors) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), layerSet.1, type=bindingFactorTypes[i], test.layer0.binding=FALSE, test.mismatch.rate=.1 ) 
  
}



head(transcriptTable)
  # one before the start co-ordinate used to get this DNA sequence.

tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.

tss.vector <- rep(0, nchar(targetSeq))
tss.vector[tss.positions] <- 1

# widen the tss positions to simulate promoters?

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



system.time(result <- optimiseFactorSet(layerList=layerList.1, factorSetRandom, testing.function=test_function, 
                                        target.layer="LAYER.1", target.vec=tss.vector, n.iter=n.iter, mut.rate=mut.rate, 
                                        modsPerCycle=modsPerCycle,logFile=logFile,logCycle=logCycle, maxNoChange=maxNoChange))

plot(result$optimScores)

save(factorSetRandom, result, file=paste("results/run.fast", n.iter,"Rdata", sep="."))

