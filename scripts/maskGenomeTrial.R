

# attempts to "mask" parts of the genome and exclude from searching for speicfic factors.
# works but does not speed up matching to genome sequence as this is always done first.

library(GenomicLayers)
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10


layerGenome <- createLayerSet.BSgenome(genome=genome, 
                                       layer.names=c("CpG_island",  "H3K27me1",  "H3K27me2","H3K27me3", "MASK"),
                                       n.layers=5)

class(layerGenome$layerSet$LAYER.0)                        
focalChrom <- "chr17"
focalLength <- seqlengths(genome)[focalChrom]
layerGenome$layerSet$MASK <- GRanges(seqnames=focalChrom, IRanges(start=1, end=focalLength), seqinfo=seqinfo(genome))

# apply a mask to layerGenome to hide all but a focal chrom   (not sure if this will work)
layerGenome

CGI.mask <- createBindingFactor.DNA_regexp("CGI", patternString="(CG.{0,4}){3}CG", patternLength=20,
                                           profile.layers = c("MASK"), profile.marks=c(1),
                                           mod.layers = "CpG_island", mod.marks=1, stateWidth=20)

bfSet <- list(CGI.mask=CGI.mask)   # binding factors must be in a list and named. Easiest to use each BFs name.

newLayerGenome <- runLayerBinding.BSgenome(layerList=layerGenome, factorSet=bfSet, iterations=100000)    # 2 mins
newLayerGenome <- runLayerBinding.BSgenome(layerList=newLayerGenome, factorSet=bfSet, iterations=100000)  # 1 second. after caching


# using a MASK layer doesn't speed up the initial search. it may save on memory.



#runLayerBinding.BSgenome(layerList = )
