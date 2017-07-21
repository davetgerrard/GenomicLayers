
# maual script to load the partial results on my PC and try out runLayerBinding and produce graphs.

# Keep running optimisation until some threshold is passed or until a set number of runs have passed without improvement
# Log results every x iterations.
# THIS VERSION:  5 layers, 1kb around the tss
# 

require(Biostrings)
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('/mnt/fls01-home01/mqbssdgb/scratch/seqPredict/predictfromsequence/scripts/pfs.functions.R')
source('scripts/pfs.functions.R') 

library(BSgenome.Hsapiens.UCSC.hg19) # note the other packages being loaded.
#available.genomes()

genome <- BSgenome.Hsapiens.UCSC.hg19
thisChrom <- genome[["chr19"]] 


# transcript.file <- '/mnt/fls01-home01/mqbssdgb/scratch/seqPredict/data/hg19.G19.chr19.transcript.gtf'
transcript.file <- "data/hg19.G19.chr19.transcript.gtf"


## Try some alternative algorithms for running a series of mods over a sequence.
# Main target is improved speed.








n.cores <- 3
base.0 <- 0 
n.layers <- 5
target.layer <- "LAYER.5"
n.factors <- 30
upstream.prom <- 200
downstream.prom <- 200
n.iter <- 10000
mut.rate <- 0.1
modsPerCycle <- 1000000 # scaled up so ~ 1/200bp
logCycle<- 100 
maxNoChange<- 1000
runName <- "pfs_layer5_chr19_1kb_par"
#outputDir <- paste("/mnt/fls01-home01/mqbssdgb/scratch/seqPredict/results/",runName, sep="")
#logFile <- paste(outputDir, "/" , runName, ".out.tab", sep="")
#outputFile <-  paste(outputDir, "/" , runName, ".final.Rdata", sep="")
#profFile <- paste(outputDir, "/" , runName, ".code.profile.out", sep="")

#if(!file.exists(outputDir)) {
#	dir.create(outputDir, recursive=TRUE)
#}

# load the (partial) results
run.result.dir <- paste("data/HYDRA_runs/", runName, "/", sep="")
final.iter <- 700
load(paste(run.result.dir, "currentFactorSet.", final.iter, ".Rdata", sep=""))

print.bfSet(currentFactorSet)


#print("creating empty layer")
#emptyLayer <-  BString(paste(rep(0,length(targetSeq)), collapse=""))


print("creating all layers")
layerSet.1 <- list(LAYER.0 = thisChrom)
for(i in 1:n.layers) {
  layerSet.1[[paste('LAYER.', i , sep="")]] <- IRanges()
}

layerList.1 <- list(layerSet=layerSet.1, history=NULL)


transcriptTable <- read.delim(transcript.file)
names(transcriptTable)[c(4,5,7)] <- c("txStart", "txEnd", "strand")

tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.

#tss.vector <- rep(0, nchar(targetSeq))
#tss.vector[tss.positions] <- 1
print(paste(sum(is.na(tss.positions)), "NAs in tss positions"))

tss.positions <- na.omit(tss.positions)
#tss.IR <- IRanges(start=tss.positions, width=1)
#tss.IR <- resize(tss.IR, fix="center", width="200")
tss.IR <- IRanges(start=tss.positions-499, end=tss.positions+500)   # this version of IRanges, resize not working right.

test_function <- function(layerList, targetLayer=target.layer, target.vec)  {
  inter.size <- sum(width(intersect(layerList$layerSet[[targetLayer]], target.vec)))
  union.size <- sum(width(union(layerList$layerSet[[targetLayer]], target.vec)))
  #layer.vec <- as.numeric(strsplit(as.character(layerList$layerSet[[targetLayer]]),"")[[1]])
  return(inter.size/ union.size)
}


# Modify the starting layer set using the optimise factorSet.

modLayerSet <- runLayerBinding(layerList=layerList.1, currentFactorSet, iterations=modsPerCycle,  collect.stats=TRUE, target.layer=target.layer, verbose=TRUE)  

modLayerSet2 <- runLayerBinding(layerList=layerList.1, currentFactorSet, iterations=modsPerCycle,  collect.stats=TRUE, target.layer=target.layer, verbose=TRUE, 
                               watch.function=function(x, target.vec=tss.IR , ...) {print(paste("Watch function", length(x[[target.layer]]), sum(overlapsAny(target.vec,x[[target.layer]]))))})  



#function(layerList, factorSet, iterations=1, bindingFactorFreqs=rep(1, length(factorSet)), watch.function=function(x){}, collect.stats=FALSE, target.layer=2, verbose=FALSE)  {

#optimiseFactorSet(layerList=layerList.1, factorSetRandom, testing.function=test_function, 
#target.layer=target.layer, target.vec=tss.IR, n.iter=n.iter, mut.rate=mut.rate, 
#modsPerCycle=modsPerCycle,logFile=logFile,logCycle=logCycle, maxNoChange=maxNoChange,
#verbose=FALSE, , use.parallel=TRUE, n.cores=n.cores))


# calc the score.
(cur.score <- test_function(modLayerSet, targetLayer=target.layer, target.vec=tss.IR))


# show the marked regions and the target regions.
# struggling with how to display whole chromosome?
# maybe something like  C:\Users\Dave\HanleyGroup\BerryChIPseq_multiTissue\scripts\regulationModelPlots.R?

#plot(start(tss.IR))

#perhaps start with a window of 1Mb

?subsetByOverlaps

start.1 <- 5000000
end.1 <- 6000000

window <- IRanges(start=start.1, end=end.1)
start(subsetByOverlaps(tss.IR, window))

set.1 <- subsetByOverlaps(tss.IR, window)
index.1 <- tileRegions(as.data.frame(set.1))
range.1 <- range(set.1)


set.2 <- subsetByOverlaps(modLayerSet$layerSet$LAYER.5,window) 
index.2 <- tileRegions(as.data.frame(set.2))
range.2 <- range(set.2)

plot(NA, xlim=c(start.1,end.1), ylim=c(0,max(index.1)))
segments(start(set.1), index.1, end(set.1), index.1, col="blue", lwd=3)
segments(start(set.1), 0, end(set.1), 0, col="black", lwd=3)

  
plot.layerAndTargets <- function(layerSet, target.layer, target.vec=NULL, win.start=0, win.end=nchar(layerSet$layerSet$LAYER.0)) {
  window <- IRanges(start=win.start, end=win.end)
  
  set.1 <- subsetByOverlaps(target.vec, window)
  index.1 <- tileRegions(as.data.frame(set.1))
  range.1 <- range(set.1)
  
  
  set.2 <- subsetByOverlaps(layerSet$layerSet[[target.layer]],window) 
  index.2 <- tileRegions(as.data.frame(set.2))
  range.2 <- range(set.2)
  
  plot(NA, xlim=c(win.start,win.end), ylim=c(0,max(index.1)))
  segments(start(set.1), index.1, end(set.1), index.1, col="blue", lwd=3)
  segments(start(set.2), 0, end(set.2), 0, col="black", lwd=3)
  
}
  
plot.layerAndTargets(layerSet=modLayerSet, target.layer=target.layer, target.vec=tss.IR, win.start=5000000, win.end=6000000)

plot.layerAndTargets(layerSet=modLayerSet, target.layer=target.layer, target.vec=tss.IR, win.start=5600000, win.end=5800000)
plot.layerAndTargets(layerSet=modLayerSet, target.layer=target.layer, target.vec=tss.IR, win.start=5750000, win.end=5800000)



# empty(ish) region
plot.layerAndTargets(layerSet=modLayerSet, target.layer=target.layer, target.vec=tss.IR, win.start=2000000, win.end=3000000)
plot.layerAndTargets(layerSet=modLayerSet, target.layer=target.layer, target.vec=tss.IR, win.start=2600000, win.end=2700000)
plot.layerAndTargets(layerSet=modLayerSet, target.layer=target.layer, target.vec=tss.IR, win.start=2620000, win.end=2660000)
#subsetByOverlaps(modLayerSet$layerSet$LAYER.5,IRanges(5780000, 5800000))


# look at the density change over the chromosome... --------------------
par(mfrow=c(2,1))
plot(density(start(tss.IR)), main="target vector")
plot(density(start(modLayerSet$layerSet$LAYER.5)), "modified target layer")
# until we have good telomere data, there may be a problem. (estimates of factor densities/occupancy)







stopifnot(FALSE)
###################

#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), n.factors, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T)
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T, prob=c(10,10,2,2))
print("generating initial factorSet")
factorSetRandom <- list()
for(i in 1:n.factors) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), layerSet.1, type=bindingFactorTypes[i], test.layer0.binding=FALSE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)

head(transcriptTable)
  # one before the start co-ordinate used to get this DNA sequence.

# widen the tss positions to simulate promoters?


test_function <- function(layerList, targetLayer=target.layer, target.vec)  {
  inter.size <- sum(width(intersect(layerList$layerSet[[targetLayer]], target.vec)))
  union.size <- sum(width(union(layerList$layerSet[[targetLayer]], target.vec)))
  #layer.vec <- as.numeric(strsplit(as.character(layerList$layerSet[[targetLayer]]),"")[[1]])
  return(inter.size/ union.size)
}



print("beginning optimisation")
#stopifnot(FALSE)

try(
system.time(result <- optimiseFactorSet(layerList=layerList.1, factorSetRandom, testing.function=test_function, 
                                        target.layer=target.layer, target.vec=tss.IR, n.iter=n.iter, mut.rate=mut.rate, 
                                        modsPerCycle=modsPerCycle,logFile=logFile,logCycle=logCycle, maxNoChange=maxNoChange,
					verbose=FALSE, , use.parallel=TRUE, n.cores=n.cores))
)
print(.Last.value)  
#plot(result$optimScores)

save(factorSetRandom, result, file= outputFile)
print("ended optimisation")
#Rprof(NULL)
#print("R code profile")
#print(summaryRprof(profFile))
#unlink(tmp)
print("end of script")

