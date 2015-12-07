

# Keep running optimisation until some threshold is passed or until a set number of runs have passed without improvement
# Log results every x iterations.
# THIS VERSION:  5 layers

require(Biostrings)
#setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')
source('/mnt/fls01-home01/mqbssdgb/scratch/seqPredict/predictfromsequence/scripts/predictFromSequence.functions.R')
source('/mnt/fls01-home01/mqbssdgb/scratch/seqPredict/predictfromsequence/scripts/pfs.functions.R')


library(BSgenome.Hsapiens.UCSC.hg19) # note the other packages being loaded.
#available.genomes()

genome <- BSgenome.Hsapiens.UCSC.hg19
thisChrom <- genome[["chr22"]] 


#fasta.file <- '/mnt/genome-shared-data/bcf/genomeIndexes/hg19_GRCh37_random_chrM/fasta/chr22.fa'
transcript.file <- '/mnt/fls01-home01/mqbssdgb/scratch/seqPredict/data/hg19.G19.chr22.transcript.gtf'

## Try some alternative algorithms for running a series of mods over a sequence.
# Main target is improved speed.

#targetSeq <- readDNAStringSet(fasta.file)[[1]]
transcriptTable <- read.delim(transcript.file)




names(transcriptTable)[c(4,5,7)] <- c("txStart", "txEnd", "strand")

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
runName <- "pfs_layer5_chr22_verbose"
outputDir <- paste("/mnt/fls01-home01/mqbssdgb/scratch/seqPredict/results/",runName, sep="")
logFile <- paste(outputDir, "/" , runName, ".out.tab", sep="")
outputFile <-  paste(outputDir, "/" , runName, ".final.Rdata", sep="")
profFile <- paste(outputDir, "/" , runName, ".code.profile.out", sep="")

if(!file.exists(outputDir)) {
	dir.create(outputDir, recursive=TRUE)
}

#print(paste(" code profiling in ", profFile))
#Rprof(profFile)    # begin profiling


print("creating empty layer")
#emptyLayer <-  BString(paste(rep(0,length(targetSeq)), collapse=""))


print("creating all layers")
layerSet.1 <- list(LAYER.0 = thisChrom)
for(i in 1:n.layers) {
  layerSet.1[[paste('LAYER.', i , sep="")]] <- IRanges()
}

layerList.1 <- list(layerSet=layerSet.1, history=NULL)



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


tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.

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



print("beginning optimisation")
#stopifnot(FALSE)

try(
system.time(result <- optimiseFactorSet(layerList=layerList.1, factorSetRandom, testing.function=test_function, 
                                        target.layer=target.layer, target.vec=tss.IR, n.iter=n.iter, mut.rate=mut.rate, 
                                        modsPerCycle=modsPerCycle,logFile=logFile,logCycle=logCycle, maxNoChange=maxNoChange,
					verbose=TRUE))
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

