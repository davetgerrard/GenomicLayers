
### simple attempt at optimising score of a layer.
require(Biostrings)
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')
source('scripts/predictFromSequence.functions.R')


#targetSeq <- DNAString('TGCGTTGAGC')
targetSeq <- readDNAStringSet("C:/Users/Dave/Dropbox/Temp/hg19.HOXA.fa")[[1]]
emptyLayer <-  BString(paste(rep(0,length(targetSeq)), collapse=""))


transcriptTable <- read.delim("data/hg19.HOXA.transcripts.tab")
head(transcriptTable)
base.0 <- 27106375   # one before the start co-ordinate used to get this DNA sequence.

tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.

tss.vector <- rep(0, nchar(targetSeq))
tss.vector[tss.positions] <- 1

n.layers <- 5

layerSet.5 <- list(LAYER.0 = targetSeq)
for(i in 1:n.layers) {
  layerSet.5[[paste('LAYER.', i , sep="")]] <- emptyLayer
}

n.factors <- 30

#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), n.factors, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T)
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T, prob=c(10,10,2,2))
factorSetRandom <- list()
for(i in 1:n.factors) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), layerSet.5, type=bindingFactorTypes[i], test.layer0.binding=FALSE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)


modLayerSet <- runLayerBinding(layerSet=layerSet.5, factorSet = factorSetRandom, iterations=100)



test_function <- function(layerSet, targetLayer="Layer.1", target.vec)  {
  layer.vec <- as.numeric(strsplit(as.character(layerSet[[targetLayer]]),"")[[1]])
  return(cor(target.vec, layer.vec))
}


test_function(modLayerSet, targetLayer="LAYER.5", target.vec=tss.vector)



# test the mutate
newFactorSet <- mutateFactorSet(factorSetRandom, layerSet.5,n.muts=3, verbose=T)


identical(newFactorSet, factorSetRandom)

print.bfSet(factorSetRandom)
print.bfSet(newFactorSet)





result <- optimiseFactorSet(layerSet.5, factorSetRandom, testing.function=test_function, target.layer="LAYER.5", target.vec=tss.vector, n.iter=500, target.score=1, mut.rate=.4)

plot(result$optimScores)

save(result, file="results/first.run.Rdata")


# TEST 2 More cycles to apply mods within each iterations
# e.g. 1000 mods, perhaps fewer iterations of mutation.


result2 <- optimiseFactorSet(layerSet.5, factorSetRandom, testing.function=test_function, target.layer="LAYER.5", target.vec=tss.vector, n.iter=100, target.score=1, mut.rate=.4, modsPerCycle=1000)

plot(result2$optimScores)

save(result2, file="results/run2.Rdata")



# RUN 3
# Have modified createRandomBindingFactor to test for binding to DNA within LAYER.0 (the DNA) of the target layerSet.
# Takes longer for initial set up and favours shorter DNA motifs. - FIXED that was due to matchPattern(,fixed=TRUE), which essentially guaranteed no matches for degenerate and/or long motifs
# should have much greater number of mods now.
# The parameter test.layer0.binding=TRUE/FALSE also propagates from optimiseFactorSet() to mutateFactorSet()

n.layers <- 1
target.layer <- "LAYER.1"

layerSet.1 <- list(LAYER.0 = targetSeq)
for(i in 1:n.layers) {
  layerSet.1[[paste('LAYER.', i , sep="")]] <- emptyLayer
}


test_function <- function(layerSet, targetLayer="Layer.1", target.vec)  {
  layer.vec <- as.numeric(strsplit(as.character(layerSet[[targetLayer]]),"")[[1]])
  return(cor(target.vec, layer.vec))
}


n.factors <- 30
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), n.factors, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T)
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T, prob=c(10,10,2,2))
factorSetRandom <- list()
for(i in 1:n.factors) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), layerSet.2, type=bindingFactorTypes[i], test.layer0.binding=TRUE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)
# picks mostly simple motifs.

modLayerSet1 <- runLayerBinding(layerSet=layerSet.1, factorSet = factorSetRandom, iterations=1000)
lapply(modLayerSet1[-1], letterFrequency, letters="1")
layer.vector1 <- as.numeric(strsplit(as.character(modLayerSet1[[target.layer]]),"")[[1]])
cor(tss.vector, layer.vector1)

gc()
# reduced the mutation rate to .1 to reduce chance of throwing baby out with bath-water.
result3 <- optimiseFactorSet(layerSet.1, factorSetRandom, testing.function=test_function, target.layer="LAYER.1", target.vec=tss.vector, n.iter=50, target.score=1, mut.rate=.1, modsPerCycle=1000, test.layer0.binding=TRUE)

plot(result3$optimScores)

save(result3, file="results/run3.Rdata")

# results very poor again.
# generate second set modLayerSet to test correlation between runs from same starting position
modLayerSet2 <- runLayerBinding(layerSet=layerSet.1, factorSet = factorSetRandom, iterations=1000)
lapply(modLayerSet2[-1], letterFrequency, letters="1")
layer.vector2 <- as.numeric(strsplit(as.character(modLayerSet2[[target.layer]]),"")[[1]])
cor(layer.vector1, layer.vector2)
# hmmm, not well correlated (perhaps should not have bothered running this). 

#  how much of the layerSet target layer is covered after each iteration of runLayerBinding?
modLayerSet3 <- runLayerBinding(layerSet=layerSet.1, factorSet = factorSetRandom, iterations=1000, 
	watch.function=function(x) print(as.numeric(unlist(lapply(x[-1], letterFrequency, letters="1")))))

#I used copy-paste to extract the coverage after each iteration (TODO add in option to return stats like this)
# It showed a monotonic increase in coverage upto 30k (of >200k), with very few instances of reduction
# Hence, with 30 such factors, and 1000 iterations, still nowhere near steady-state or saturation.
# Need ~10 times as many iterations 
#  (and ability to store 
# TODO re-design layerSet to have list of layers at lower level and have a top level with optimisation records.
#		e.g. n.iterations, target.score, layer.states per iteration. Perhaps even the bindingFactor name?



n.layers <- 1
target.layer <- "LAYER.1"

layerSet.1 <- list(LAYER.0 = targetSeq)
for(i in 1:n.layers) {
  layerSet.1[[paste('LAYER.', i , sep="")]] <- emptyLayer
}

layerList.1 <- list(layerSet=layerSet.1, history=NULL)

layerList.1$layerSet
# THIS MODIFICATION BREAKS EVERYTHING THAT HAS GONE BEFORE...
modLayerSet4 <- runLayerBinding(layerList=layerList.1, factorSet = factorSetRandom, iterations=10000, collect.stats=TRUE)
plot(modLayerSet4$history$target.coverage)  # after 10000 mods this 200kb region with 1-layer was beginning to level off. Going to need a huge number of mods genome wide!

save(modLayerSet4 , file="results/modLayerSet4.10000.stats.Rdata")

table(modLayerSet4$history$bf)  # most of the bindingFactors got used an equal number of times. Except for large layer_island types (e.g. 0000111110000) which will have struggled to find a match.


# TODO
# Test with different sequence.
# Test using a simpler base sequence with known motifs strongly tied to the output.  i.e. how hard is to find a motif when there is one?
# Test using the same starter sequence but with biologically inspired motifs e.g. TATA, GC-rich.
# Enable recording of a sequence of binding events so that it can be 're-layed' as a graphic and stats collected about the modifications applied by each factor (and which factors have no effect)
# Move testing to Kadmon or Hydra with email notification on finish
# Allow for fixed modificatin order?  No, I want to emulate the statistical randomness of factor proportions.
# Allow for varying proportions of components? 
# Pre-generate list of random 'matching' bindingFactors on whole genome. (current strategy to test 1000 is failing to find a match for DNA_region types and for DNA_motifs > 8.
# N.B. prediction of nucleosome positions:- http://genie.weizmann.ac.il/software/nucleo_exe.html



