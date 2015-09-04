
require(Biostrings)
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')
source('scripts/predictFromSequence.functions.R')

# N.B. orientation of transcripts may be a problem

# load HOXD data

#targetSeq <- readDNAStringSet("C:/Users/Dave/Dropbox/Temp/hg19.HOXA.fa")[[1]]
targetSeq <- readDNAStringSet("data/hg19.HOXD.fa")[[1]]
emptyLayer <-  BString(paste(rep(0,length(targetSeq)), collapse=""))



n.layers <- 5

layerSet <- list(LAYER.0 = targetSeq)
for(i in 1:n.layers) {
  layerSet[[paste('LAYER.', i , sep="")]] <- emptyLayer
}

layerList.5 <- list(layerSet=layerSet, history=NULL)

# load transcript table and create tss.vector

#transcriptTable <- read.delim("data/hg19.HOXA.transcripts.tab")
transcriptTable <- read.delim("data/hg19.HOXD.transcripts.tab")
head(transcriptTable)
base.0 <- 176893696   # one before the start co-ordinate used to get this DNA sequence.

tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.
length(tss.positions)  # unique tss
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


# load factorSet optimised on HOXA
load("data//HYDRA_runs//layer5_10k//layer5_10k.final.Rdata")

targetLayer <- "LAYER.1"

# run the initial factor set mods 10 times to get distribution
#initialMod <- 
# use more, don't keep the full result calc sum along the way.
reps <- 100
result.vector.initial <- rep(0, length(tss.vector))
for(i in 1:reps) {
  print(".")
  thisModLayerList <- runLayerBinding.fast(layerList=layerList.5, factorSet = factorSetRandom, iterations=10000)  # result[-(length(result))] because last object is a set of scores
  this.result.vector <- as.numeric(strsplit(as.character(thisModLayerList$layerSet[[targetLayer]]),"")[[1]])
  result.vector.initial <- result.vector.initial + this.result.vector
}

# repeat with final factor set to get distribution
result.vector.final <- rep(0, length(tss.vector))
for(i in 1:reps) {
  print(".")
  thisModLayerList <- runLayerBinding.fast(layerList=layerList.5, factorSet = result[1:30], iterations=10000)  # result[-(length(result))] because last object is a set of scores
  this.result.vector <- as.numeric(strsplit(as.character(thisModLayerList$layerSet[[targetLayer]]),"")[[1]])
  result.vector.final <- result.vector.final + this.result.vector
}  
  
  
# plot them 
plot(which(tss.vector == 1), rep(1.05, sum(tss.vector)), xlim=c(0, length(tss.vector)), ylim=c(0,1.1), xlab="position", ylab="Proportion of mods base is altered")
points(result.vector.initial/reps, col="blue")
points(which(tss.vector == 1), (result.vector.initial/reps)[tss.vector == 1], col="green")

plot(which(tss.vector == 1), rep(1.05, sum(tss.vector)), xlim=c(0, length(tss.vector)), ylim=c(0,1.1), xlab="position", ylab="Proportion of mods base is altered")
points(result.vector.final/reps, col="blue")
points(which(tss.vector == 1), (result.vector.final/reps)[tss.vector == 1], col="green")

# scores
cor(tss.vector, result.vector.initial/reps)
cor(tss.vector, result.vector.final/reps)
# any improvement?   - yes
  

save(result.vector.initial, result.vector.final, file="results/layer_5_test.100modsSummed.Rdata")
  
  