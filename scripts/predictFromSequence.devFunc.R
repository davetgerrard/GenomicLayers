## Scratchbox to test the latest version of the functions.


# Investigating low coverage on optimisation --------------------
# Increase coverage of modifications?
require(Biostrings)
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
source('scripts/pfs.functions.R')   


# What to do about overlapping targets (after expansion)?

# How to increase coverage? 

# Could allow optimisation to choose coverage based on bindingFactor lengths

# Check that they are binding as much as they could... e.g.
#   bf.2 6878 1bf.3 1 1[1] "Sequence of length  51304566 , using  6 windows of length 1e+07"
# bf.4 17615 1[1] "Sequence of length  51304566 , using  6 windows of length 1e+07"
# bf.7 37536 1[1] "Sequence of length  51304566 , using  6 windows of length 1e+07"

library(BSgenome.Hsapiens.UCSC.hg19) # note the other packages being loaded.
#available.genomes()

genome <- BSgenome.Hsapiens.UCSC.hg19
thisChrom <- genome[["chrM"]] 
n.layers <- 5

layerSet.5 <- list(LAYER.0 = thisChrom)
for(i in 1:n.layers) {
  layerSet.5[[paste('LAYER.', i , sep="")]] <- IRanges()    # use IRanges to store state of layers. TODO limit to chrom length
}

layerList.5 <- list(layerSet=layerSet.5, history=NULL)

# generate a factor set and test how many alter the target layer and by how much.

n.factors <- 30
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), n.factors, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T)
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T, prob=c(10,10,2,2))
factorSetRandom <- list()
for(i in 1:n.factors) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), 
                                                                         layerSet.5, type=bindingFactorTypes[i], test.layer0.binding=TRUE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)

# create a simple binding factor that binds anywhere and convert the target layer (LAYER.5) to state 1
testFactor <-  list(name="test", type="layer_region", 
                                   profile=list(LAYER.0=list(pattern=DNAString("N"), mismatch.rate=0.1, length=250)), 
                                   mods=list(LAYER.5=list(state="1", stateWidth=250, offset=0, align="centre")))

# get hits for a specific factor

hits <- matchBindingFactor(layerSet=layerSet.5, bindingFactor=testFactor, verbose=T)   ## TODO - get this working. 
#hits <- as(hits, "IRanges")  # needed to use reduce, setDiff etc TODO fix matchBindingFactor to return this?
print( paste(testFactor$type,length(hits), class(hits), "hits"))
#print(hits)
new.LayerSet <- modifyLayerByBindingFactor.Views(layerSet.5, hits=hits, bindingFactor=testFactor)


modLayerSet <- runLayerBinding(layerList=layerList.5, factorSet = list(A=testFactor), iterations = 10000, verbose=TRUE)

# attempt optimisation using parallel runs on different mutant sets --------------------------
require(Biostrings)
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')
#source('scripts/predictFromSequence.functions.R')
source('scripts/pfs.functions.R')   # overwrites some of the above. TODO - remove this dependency.

library(BSgenome.Hsapiens.UCSC.hg19) # note the other packages being loaded.
#available.genomes()

genome <- BSgenome.Hsapiens.UCSC.hg19
thisChrom <- genome[["chrM"]] 
transcript.file <- 'data/hg19.G19.chrM.transcript.gtf'
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
runName <- "pfs_layer5_chr22_1kb_verbose"



layerSet.5 <- list(LAYER.0 = thisChrom)
for(i in 1:n.layers) {
  layerSet.5[[paste('LAYER.', i , sep="")]] <- IRanges()    # use IRanges to store state of layers. TODO limit to chrom length
}

layerList.5 <- list(layerSet=layerSet.5, history=NULL)

n.factors <- 30
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), n.factors, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T)
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T, prob=c(10,10,2,2))
factorSetRandom <- list()
for(i in 1:n.factors) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), 
                                                                         layerSet.5, type=bindingFactorTypes[i], test.layer0.binding=TRUE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)

## Try some alternative algorithms for running a series of mods over a sequence.
# Main target is improved speed.

#targetSeq <- readDNAStringSet(fasta.file)[[1]]
transcriptTable <- read.delim(transcript.file)
names(transcriptTable)[c(4,5,7)] <- c("txStart", "txEnd", "strand")

# set up for optimisation
# target.vec is now an IRanges
tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions)    # only want unique values.

#tss.vector <- rep(0, nchar(targetSeq))
#tss.vector[tss.positions] <- 1
print(paste(sum(is.na(tss.positions)), "NAs in tss positions"))

tss.positions <- na.omit(tss.positions)
#tss.IR <- IRanges(start=tss.positions, width=1)
#tss.IR <- resize(tss.IR, fix="center", width="200")
tss.IR <- IRanges(start=tss.positions-99, end=tss.positions+100)   # this version of IRanges, resize not working right.

# widen the tss positions to simulate promoters?


test_function <- function(layerList, targetLayer=target.layer, target.vec)  {
  inter.size <- sum(width(intersect(layerList$layerSet[[targetLayer]], target.vec)))
  union.size <- sum(width(union(layerList$layerSet[[targetLayer]], target.vec)))
  #layer.vec <- as.numeric(strsplit(as.character(layerList$layerSet[[targetLayer]]),"")[[1]])
  return(inter.size/ union.size)
}

# TODO : get these object back into aligment.

# test the Layer binding
system.time(modLayerSet.fast <- runLayerBinding(layerList=layerList.5, factorSet = factorSetRandom, verbose=TRUE))  


n.iter <- 100
#system.time(result <- optimiseFactorSet(layerList=layerList.1, factorSetRandom, testing.function=test_function, target.layer="LAYER.5", target.vec=tss.IR, n.iter=n.iter, mut.rate=0.1, modsPerCycle=10000, verbose=TRUE))

try(
  system.time(result <- optimiseFactorSet(layerList=layerList.5, factorSetRandom, testing.function=test_function,
                                          target.layer=target.layer, target.vec=tss.IR, n.iter=n.iter, mut.rate=mut.rate,
                                          modsPerCycle=modsPerCycle,logFile="temp.log",logCycle=logCycle, maxNoChange=maxNoChange,
                                          verbose=TRUE))
)


#mutateAndLayerBindFactorSet <- function(factorSet, layerList, mut.rate =0.1,  verbose=TRUE, test.layer0.binding=FALSE, name.prefix="m", modsPercyle=100)  {
#  n.muts=floor(length(currentFactorSet) * mut.rate)
#  
#  newFactorSet <- mutateFactorSet(factorSet, layerList$layerSet,n.muts=n.muts, verbose=verbose, test.layer0.binding=test.layer0.binding, name.prefix=name.prefix)
#    
#  newModLayer <- runLayerBinding(layerList=layerList, factorSet = newFactorSet, iterations=modsPerCycle, verbose=verbose)
#  return(newModLayer)#
#}

#mutateAndLayerBindFactorSet(factorSet, layerList, mut.rate =0.1,  verbose=TRUE, test.layer0.binding=TRUE, name.prefix=paste("m",i,sep=""), modsPercyle=100)

#newFactorSet <- mutateFactorSet(currentFactorSet, layerList$layerSet,n.muts=floor(length(currentFactorSet) * mut.rate), verbose=T, test.layer0.binding=test.layer0.binding, name.prefix=paste("m",i,sep=""))
#newModLayer <- runLayerBinding(layerList=layerList, factorSet = newFactorSet, iterations=modsPerCycle, verbose=verbose)

# NO, need it to return the newModLayer AND keep track of which factorSet went in.
# So, create the N factor sets
# then do runLayerBinding in parallel

factorSet.list <- list()

#mutateFactorSet(factorSetRandom, layerList.5$layerSet,n.muts=floor(length(factorSetRandom) * mut.rate), verbose=T, test.layer0.binding=T, name.prefix=paste("m",i,sep=""))


# create 10 new factorSet by mutating the currentSet.
factorSet.list <- replicate(10, mutateFactorSet(factorSetRandom, layerList.5$layerSet,n.muts=floor(length(factorSetRandom) * mut.rate), verbose=T, test.layer0.binding=T, name.prefix=paste("m",i,sep="")), simplify=FALSE)

#length(factorSet.list)
#length(factorSet.list[[1]])
#print.bfSet(factorSet.list[[1]])

#sapply(factorSet.list, print.bfSet, simplify=F)

# attempt at parallel runs
require(parallel)

cl <- makeCluster(mc <- getOption("cl.cores", 3))
## to make this reproducible
#clusterSetRNGStream(cl, 123)

clusterExport(cl, c("runLayerBinding", "matchBindingFactor","modifyLayerByBindingFactor.Views", "layerList.5", "modsPerCycle", "verbose"))
n.tries <- 10  # VERY WEIRDLY , the below command sometimes fails on first call but works on second (or later).  TODO: fix this!
tries=0
while(tries <= n.tries)  {
  
  worked <- try(mod.list <- parLapply(cl, factorSet.list, fun=function(x) runLayerBinding(layerList=layerList.5, factorSet = x, iterations=modsPerCycle, verbose=verbose) ), TRUE)
  tries <- tries +1
  if(class(worked) != "try-error") {
    if(verbose) print(paste("parallel mod list tries: " , tries))
    break;
  }
}
clusterExport(cl, c("test_function", "target.layer","tss.IR"))
score.vec <- parSapply(cl, mod.list, FUN=function(x) test_function(layerList=x, targetLayer=target.layer, target.vec=tss.IR))
stopCluster(cl)




# TESTING
# simple 60 iters
try(
  system.time(result <- optimiseFactorSet(layerList=layerList.5, factorSetRandom, testing.function=test_function,
                                          target.layer=target.layer, target.vec=tss.IR, n.iter=60, mut.rate=mut.rate,
                                          modsPerCycle=modsPerCycle,logFile="temp.log",logCycle=logCycle, maxNoChange=maxNoChange,
                                          verbose=TRUE))
)
# RUN1
# [1] "Round 60 oldScore 0.0622369878183832 newScore 0.0321151716500554 Better?"
# [1] "No!"
# user  system elapsed 
# 94.13    0.07   95.19 
# RUN2
#[1] "Round 60 oldScore 0.1828125 newScore 0.154369440083726 Better?"
#[1] "No!"
#user  system elapsed 
#112.56    0.02  113.17 

# source('scripts/pfs.functions.R')  
# 20 iters using 3 cores
try(
  system.time(result <- optimiseFactorSet(layerList=layerList.5, factorSetRandom, testing.function=test_function,
                                          target.layer=target.layer, target.vec=tss.IR, n.iter=20, mut.rate=mut.rate,
                                          modsPerCycle=modsPerCycle,logFile="temp.log",logCycle=logCycle, maxNoChange=maxNoChange,
                                          verbose=T, use.parallel=TRUE, n.cores=3))
)
# RUN1
#[1] "Round 20 oldScore 0.175549322113137 newScore 0.175722947761194 Better?"
#[1] "Yes!"
#user  system elapsed 
#3.71    0.46   52.28 
# RUN2
# [1] "Round 20 oldScore 0.184375399565273 newScore 0.175470239628962 Better?"
# [1] "No!"
# user  system elapsed 
# 4.26    0.39   61.65 
# RUN3
# [1] "Round 20 oldScore 0.0680616740088106 newScore 0.0680616740088106 Better?"
# [1] "No!"
# user  system elapsed 
# 3.84    0.49   54.04

# Three cores on my PC is faster and at least as good.

# overnight run, just to see how far it will go on chrM
try(
  system.time(result <- optimiseFactorSet(layerList=layerList.5, factorSetRandom, testing.function=test_function,
                                          target.layer=target.layer, target.vec=tss.IR, n.iter=10000, mut.rate=mut.rate,
                                          modsPerCycle=modsPerCycle,logFile="temp.log",logCycle=logCycle, maxNoChange=maxNoChange,
                                          verbose=F, use.parallel=TRUE, n.cores=3))
)




# Test binding of multi-layers ---------------------------
require(Biostrings)
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')
source('scripts/predictFromSequence.functions.R')
source('scripts/pfs.functions.R')   # overwrites some of the above. TODO - remove this dependency.

library(BSgenome.Hsapiens.UCSC.hg19) # note the other packages being loaded.
#available.genomes()

genome <- BSgenome.Hsapiens.UCSC.hg19
thisChrom <- genome[["chrM"]] 
n.layers <- 5

layerSet.5 <- list(LAYER.0 = thisChrom)
for(i in 1:n.layers) {
  layerSet.5[[paste('LAYER.', i , sep="")]] <- IRanges()    # use IRanges to store state of layers. TODO limit to chrom length
}

layerList.5 <- list(layerSet=layerSet.5, history=NULL)

n.factors <- 30
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), n.factors, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T)
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T, prob=c(10,10,2,2))
factorSetRandom <- list()
for(i in 1:n.factors) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), 
                                                                         layerSet.5, type=bindingFactorTypes[i], test.layer0.binding=TRUE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)

new.LayerSet <- layerSet.5
for(i in 1:length(factorSetRandom)) {
  print(i)
  hits <- matchBindingFactor(layerSet=new.LayerSet, bindingFactor=factorSetRandom[[i]], verbose=T)   ## TODO - get this working. 
  #hits <- as(hits, "IRanges")  # needed to use reduce, setDiff etc TODO fix matchBindingFactor to return this?
  print( paste(factorSetRandom[[i]]$type,length(hits), class(hits), "hits"))
  #print(hits)
  new.LayerSet <- modifyLayerByBindingFactor.Views(new.LayerSet, hits=hits, bindingFactor=factorSetRandom[[i]])
  #print(new.LayerSet)
}

#  won't find very matches with 5 layers very easily.



# complete re-write to use Views or Iranges as layer matches (e.g. on BSgenome)   ---------------------------------------------
# functions in pfs.functions.R
require(Biostrings)
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')
source('scripts/predictFromSequence.functions.R')
source('scripts/pfs.functions.R')   # overwrites some of the above. TODO - remove this dependency.

library(BSgenome.Hsapiens.UCSC.hg19) # note the other packages being loaded.
#available.genomes()

genome <- BSgenome.Hsapiens.UCSC.hg19
thisChrom <- genome[["chrM"]] 
#thisChrom <- genome[["chr22"]]   # BEWARE, memory usage.

z <- Views(thisChrom)

n.layers <- 1

layerSet.1 <- list(LAYER.0 = thisChrom)
for(i in 1:n.layers) {
  layerSet.1[[paste('LAYER.', i , sep="")]] <- IRanges()    # use IRanges to store state of layers. TODO limit to chrom length
}

layerList.1 <- list(layerSet=layerSet.1, history=NULL)

n.factors <- 30
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region"), n.factors, replace=T) 
bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T)
#bindingFactorTypes <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"), n.factors, replace=T, prob=c(10,10,2,2))
factorSetRandom <- list()
for(i in 1:n.factors) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), 
                            layerSet.1, type=bindingFactorTypes[i], test.layer0.binding=TRUE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)

# this is where MEMORY problems begin
#hits <- system.time(matchBindingFactor(layerSet=layerSet.1, bindingFactor=factorSetRandom[[4]])   ## TODO - get this working. )
#hits <- as(hits, "IRanges")  # needed to use reduce, setDiff etc TODO fix matchBindingFactor to return this?
rm(hits)
system.time(hits <- matchBindingFactor(layerSet=layerSet.1, bindingFactor=factorSetRandom[[4]], max.window=10000000, verbose=T))

tail(sort(width(hits))) # using max.window seems to produce a small number of very long hits. Weird.  Could remove them with a width filter but nnot ure how they are made

new.LayerSet <- modifyLayerByBindingFactor.Views(layerSet.1, hits=hits, bindingFactor=factorSetRandom[[4]])

new.LayerSet <- layerSet.1
for(i in 1:length(factorSetRandom)) {
  print(i)
  hits <- matchBindingFactor(layerSet=new.LayerSet, bindingFactor=factorSetRandom[[i]], verbose=T)   ## TODO - get this working. 
  #hits <- as(hits, "IRanges")  # needed to use reduce, setDiff etc TODO fix matchBindingFactor to return this?
  print( paste(factorSetRandom[[i]]$type,length(hits), class(hits), "hits"))
  #print(hits)
  new.LayerSet <- modifyLayerByBindingFactor.Views(new.LayerSet, hits=hits, bindingFactor=factorSetRandom[[i]])
  print(new.LayerSet)
}



sum(width(new.LayerSet$LAYER.1))
system.time(modLayerSet.fast <- runLayerBinding(layerList=layerList.1, factorSet = factorSetRandom, verbose=TRUE))  
sum(width(modLayerSet.fast$LAYER.1))
#layerList.m <- list(layerSet=modLayerSet.fast, history=NULL)


# set up for optimisation
# target.vec is now an IRanges
#target.IR <- IRanges(start=c(3000, 5500), end=c(3500, 6000))
target.IR <- new.LayerSet$LAYER.1[sample(1:length(new.LayerSet$LAYER.1), 2000)]
test_function <- function(layerList, targetLayer="Layer.1", target.vec)  {
  inter.size <- sum(width(intersect(layerList$layerSet[[targetLayer]], target.vec)))
  union.size <- sum(width(union(layerList$layerSet[[targetLayer]], target.vec)))
  #layer.vec <- as.numeric(strsplit(as.character(layerList$layerSet[[targetLayer]]),"")[[1]])
  return(inter.size/ union.size)
}
test_function(layerList=modLayerSet.fast$layerSet, targetLayer="LAYER.1", target.vec=target.IR)

# TODO : get these object back into aligment.


n.iter <- 10
system.time(result <- optimiseFactorSet(layerList=layerList.1, factorSetRandom, testing.function=test_function, target.layer="LAYER.1", target.vec=target.IR, n.iter=n.iter, mut.rate=0.1, modsPerCycle=10000, verbose=TRUE))

modLayerSet.result <- runLayerBinding(layerList=layerList.1, factorSet = result[1:30])
test_function(layerList=modLayerSet.result, targetLayer="LAYER.1", target.vec=target.IR)


# improving speed in modifyLayers... --------------------------

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


#system.time(modLayerSet <- runLayerBinding.fast(layerList=layerList.1, factorSet = factorSetRandom))   # 10secs

x <- sample(1:10000, 50)

y <- pmin(x + 500, 10000)   # truncate to max value
y

#ll <- Views(layerSet.1$LAYER.1, start=x, end=y)
# need to create set of ranges 
# then reduce() to get non-overlapping range
# use this to change state of layer.

rr <- IRanges(start=x, end=y)
sum(width(rr))

redr <- reduce(rr)
sum(width(redr))

layerSet.1$LAYER.1[redr] <- BString('1')



# having problems with chromosome level mapping
# try BSgenome
library(BSgenome.Hsapiens.UCSC.hg19) # note the other packages being loaded.
available.genomes()

genome <- BSgenome.Hsapiens.UCSC.hg19
thisChrom <- genome[["chr22"]] 
pattern <- DNAString("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")

hits <- matchPattern(pattern, thisChrom, fixed=F) # a few seconds. 


targetSeq <- thisChrom
layerSet.x <- list(LAYER.0 = targetSeq)
emptyLayer <-  BString(paste(rep(0,length(targetSeq)), collapse=""))  # this is slow for whole chrom. Uses quite a bit of memory.
for(i in 1) {
  layerSet.x[[paste('LAYER.', i , sep="")]] <- emptyLayer
}


thisBf <- createRandomBindingFactor(layerSet=layerSet.x, name="test", type="DNA_region")
thisBf$profile$LAYER.0$pattern <- pattern
# the following line is probably why scripts on chr7 and chr22 were crashing out so hard. 
#hits.2 <- matchBindingFactor(layerSet=layerSet.x, bindingFactor=thisBf, clusterGap=10)   # MEMORY MONSTER!

# TODO: use BSgenome and views for pattern matching on the sequence
# TODO: find an alternative to BString for pattern matching and modification on the layers.
#         BSgenome?  rle? or views/hits 
#     Views/hits may work if just record regions that are marked or not. 
#           Except when trying to match islands (but how common is that anyway?) and could be done with setdiff() or gaps().
#           views/hits would be very easy to modify.

# from ?IRanges
IRanges()  # IRanges instance of length zero
## With logical input:
x <- IRanges(c(FALSE, TRUE, TRUE, FALSE, TRUE))  # logical vector input
isNormal(x)  # TRUE
x <- IRanges(Rle(1:30) %% 5 <= 2)  # logical Rle input   # how does that work?
isNormal(x)  # TRUE


# to match a certain width on IRanges
# for match to postiive (1) values
x <- IRanges(start=c(100, 200, 500), width=c(75,10,50))
x[width(x) > 20]
gaps(x)
# for match to negative (0) values
gaps(x)[width(gaps(x)) > 50]    # really need to specify upper end of range
x

# to match an island...
# one positive match adjacent to two gap matches
pos <- x[width(x) > 8]
neg <- gaps(x)[width(gaps(x)) > 20] 
countOverlaps(pos,neg, maxgap = 1)   # values of 2 are positive regions overlapping two negative regions by 1 bp (i.e. on either side).
# would need to be inverted for negative islands.



# what was this? -----------------------------------------------------


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
table.factorSet(currentFactorSet)
type.patterns(table.factorSet(currentFactorSet))


# plot a euler grid of the mod and target combinations to see what they look like.
source("c:/users/dave/utilsGerrardDT/dataToEulerGrid.R")
type.pattern.table <- as.data.frame(type.patterns(table.factorSet(currentFactorSet)))
type.pattern.table <- type.pattern.table[, sort(names(type.pattern.table))]
euler.table <- scoreCardinalities(type.pattern.table)
png("figures/test.euler.layer5_partial.png")
par(mar=c(3,10,3,3))
plotEuler(subset(euler.table, select=- count), counts=euler.table$count)
dev.off()
