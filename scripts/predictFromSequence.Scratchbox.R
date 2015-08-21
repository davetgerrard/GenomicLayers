#   Make multiple random bindingFactor sets and test each one.
require(Biostrings)
source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')



#targetSeq <- DNAString('TGCGTTGAGC')
targetSeq <- readDNAStringSet("C:/Users/Dave/Dropbox/Temp/hg19.HOXA.fa")[[1]]
emptyLayer <-  BString(paste(rep(0,length(targetSeq)), collapse=""))

n.layers <- 5

layerSet.5 <- list(LAYER.0 = targetSeq)
for(i in 1:n.layers) {
  layerSet.5[[paste('LAYER.', i , sep="")]] <- emptyLayer
}

bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), 50, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), 50, replace=T)
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), 50, replace=T, prob=c(10,10,2,2))
factorSetRandom <- list()
for(i in 1:50) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), layerSet.5, type=bindingFactorTypes[i], test.layer0.binding=FALSE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)


modLayerSet <- runLayerBinding(layerSet=layerSet.5, factorSet = factorSetRandom, iterations=100)
#modLayerSet <- runLayerBinding(layerSet=layerSet.5, factorSet = factorSetRandom, iterations=100)
gc()  # seems to use up a tonne of memory..
lapply(modLayerSet[-1], letterFrequency, letters="1")
# Was too slow when matching very long N strings for layer_region and layer_island types. Solved by reducing to a single N


# seems very slow now I've added in larger region match types. 
# suspect the search space is larger
# or perhaps the overlaps are slower to compute, or larger.

n.runs <- 100
layerSet.list <- list()
for(i in 1:n.runs) {
  layerSet.list[[i]] <- layerSet.5
}

library(parallel)

cl <- makeCluster(getOption("cl.cores", 3))
clusterExport(cl, c("runLayerBinding","modifyLayerByBindingFactor", "modifyLayer", "truncateValueWithinRange", "matchBindingFactor" ))
clusterExport(cl, c("layerSet.5", "factorSetRandom"))
#result <- clusterCall(cl, get("runLayerBinding"), layerSet=layerSet.5, factorSet = factorSetRandom, iterations=1000)
#result <- clusterApply(cl, x=1:10,  get("runLayerBinding"), layerSet=layerSet.5, factorSet = factorSetRandom, iterations=100)
result <- parLapply(cl, layerSet.list, fun = function(x) runLayerBinding(x,  factorSet = factorSetRandom, iterations=100))
#result <- parLapply(cl, layerSet.list, fun = function(x) runLayerBinding(x,  factorSet = factorSetRandom, iterations=1000))
stopCluster(cl)



letterFrequency(result[[1]]$LAYER.1, letters= "1")

lapply(result[[1]][-1], letterFrequency, letters="1")
lapply(result[[2]][-1], letterFrequency, letters="1")
lapply(result[[3]][-1], letterFrequency, letters="1")

transcriptTable <- read.delim("C:/Users/Dave/Dropbox/Temp/hg19.HOXA.transcripts.tab")
head(transcriptTable)
base.0 <- 27106375   # one before the start co-ordinate used to get this DNA sequence.

tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.

tss.vector <- rep(0, nchar(targetSeq))
tss.vector[tss.positions] <- 1



targetLayer <- 'LAYER.5'
par(mfrow=c(3,1))
cor.vec <- numeric()
for(i in 1:length(result)) {
  layer.vector <- as.numeric(strsplit(as.character(result[[i]][[targetLayer]]),"")[[1]])
  #plot(layer.vector, type="l")  
  cor.vec[i] <- cor(tss.vector, layer.vector)
}
cor.vec

hist(cor.vec)
# the results can, currently be quite variable (with 100 iterations) as the resulting marks are very sparse.

#plot(tss.vector, xlim=c(0, nchar(targetSeq)), type="l")

for(i in 1:length(result)) {
  layer.vector <- as.numeric(strsplit(as.character(result[[i]][[targetLayer]]),"")[[1]])
  plot(layer.vector, type="l")   
}



####



# try to mutate a factorSet to improve the match

# set base0 factorSet

# runLayerBnding

# get score and set as base0

# mutate factorSet
# runLayerBinding
# get score
# if score improved, reset base0 and repeat.


mutate.DNA_motif <- function(x)  {
  
}


mutate.DNA_region <- function(x)  {
  
}

mutate.layer_region <- function(x)  {
  
}

mutate.layer_island <- function(x)  {
  
}


mutateLayerMod <- function(x)  {
  
}



# master function to control mutations across a factorSet
# Mutations options:-
#   1. mutate individual bindingFactors
#   2. replace binding factorswith a new random bindingFactor
#   3. change the proportions of existing bindingFactors.  (need to keep a vector of bindingFactor frequncies)
mutateFactorSet <- function(x)  {
  
}



# need to trial evolving with much smaller datasets.




#how many of the random patterns actually match the dna?

pattern.vec <- character()
match.lengths <- integer()
for(i in 1:100) {
  patternLength <- max(1,rnbinom(1, 10, mu= 10))  # patterns must be at least of length 1
  pattern.vec[i] <- paste(sample(names(IUPAC_CODE_MAP), patternLength, replace=T), collapse="")
  match.lengths[i] <- length( matchPattern(pattern.vec[i], targetSeq))
  
}

# Wow, hardly any!  
pattern.vec[match.lengths > 0]


# now using modified function testing for matches

pattern.vec <- character()
match.lengths <- integer()
for(i in 1:30) {
  newFactor <- createRandomBindingFactor(paste("bf.",i ,sep=""), layerSet.5, type=bindingFactorTypes[i], test.layer0.binding=TRUE, test.mismatch.rate=.1 ) 
  match.lengths[i] <- length( matchBindingFactor(layerSet.5, newFactor))
  
}



############################# DEVELOPMENT
?"parallel"
#
library(parallel)
## Use option cl.cores to choose an appropriate cluster size.
cl <- makeCluster(getOption("cl.cores", 2))



clusterApply(cl, 1:2, get("+"), 3)
xx <- 1
clusterExport(cl, "xx")
result <- clusterCall(cl, function(y) xx + y, 2)
stopCluster(cl)



testFunction <- function(xx, y=3, z) {
   return(xx + y + z) 
}

cl <- makeCluster(getOption("cl.cores", 2))
result <- clusterApply(cl, x=1:10,  get("testFunction"), z=1)
stopCluster(cl)



