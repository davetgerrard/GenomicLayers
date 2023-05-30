

# method to subset a BSgenome for testing.   Speeds up searches considerably. 

library(GenomicLayers)
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10

# hack from https://support.bioconductor.org/p/83588/  to keep only some chromosomes within a BSgenome object.
keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

#sequences_to_keep <- paste0("chr", c(1:20, "X", "Y"))
sequences_to_keep <- "chr17"
genomeSub <- keepBSgenomeSequences(genome, sequences_to_keep)
genomeSub    # this should now still be a useable BSgenome object but with only one chromosome.  

CGI<- createBindingFactor.DNA_regexp("CGI", patternString="(CG.{0,4}){3}CG", patternLength=20,
                                   mod.layers = "CpG_island", mod.marks=1, stateWidth=20)

bfSet <- list(CGI=CGI)   # binding factors must be in a list and named. Easiest to use each BFs name.


layerFullGenome <- createLayerSet.BSgenome(genome=genome, 
                                       layer.names=c("CpG_island",  "H3K27me1",  "H3K27me2","H3K27me3"),
                                       n.layers=4)
layerPartialGenome <- createLayerSet.BSgenome(genome=genomeSub, 
                                           layer.names=c("CpG_island",  "H3K27me1",  "H3K27me2","H3K27me3"),
                                           n.layers=4)

# can use system.time() to measure the speed difference of one chromosome vs the whole genome
system.time(
  newLayerFullGenome <- runLayerBinding.BSgenome(layerList=layerFullGenome, factorSet=bfSet, iterations=100000) 
)
system.time(
  newLayerPartialGenome <- runLayerBinding.BSgenome(layerList=layerPartialGenome, factorSet=bfSet, iterations=100000) 
)
