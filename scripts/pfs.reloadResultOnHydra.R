
# re-load a saved run and run runLayerBinding to get some stats

library(Biostrings)
source("predictfromsequence/scripts/pfs.functions.R")
load("results/pfs_layer5_chr22_1kb_par/currentFactorSet.600.Rdata")

library(BSgenome.Hsapiens.UCSC.hg19) # note the other packages being loaded.
genome <- BSgenome.Hsapiens.UCSC.hg19
thisChrom <- genome[["chr22"]]

transcript.file <- '/mnt/fls01-home01/mqbssdgb/scratch/seqPredict/data/hg19.G19.chr22.transcript.gtf'

## Try some alternative algorithms for running a series of mods over a sequence.
# Main target is improved speed.

#targetSeq <- readDNAStringSet(fasta.file)[[1]]
transcriptTable <- read.delim(transcript.file)

names(transcriptTable)[c(4,5,7)] <- c("txStart", "txEnd", "strand")

base.0 <- 0
n.layers <- 5
target.layer <- "LAYER.5"

n.iter <- 10000
mut.rate <- 0.1
modsPerCycle <- 1000000 # scaled up so ~ 1/200bp
logCycle<- 100
maxNoChange<- 1000


print("creating all layers")
layerSet.1 <- list(LAYER.0 = thisChrom)
for(i in 1:n.layers) {
  layerSet.1[[paste('LAYER.', i , sep="")]] <- IRanges()
}

layerList.1 <- list(layerSet=layerSet.1, history=NULL)

tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.

#tss.vector <- rep(0, nchar(targetSeq))
#tss.vector[tss.positions] <- 1
print(paste(sum(is.na(tss.positions)), "NAs in tss positions"))

tss.positions <- na.omit(tss.positions)
#tss.IR <- IRanges(start=tss.positions, width=1)
#tss.IR <- resize(tss.IR, fix="center", width="200")
tss.IR <- IRanges(start=tss.positions-499, end=tss.positions+500)   # this version of IRanges, resize not working right.

modLayerSet <- runLayerBinding(layerList=layerList.1, factorSet = currentFactorSet, verbose=TRUE)


countOverlaps(modLayerSet$layerSet$LAYER.5, tss.IR)

# the answer was pretty awful.
