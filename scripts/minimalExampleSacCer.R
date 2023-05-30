# minimal example on whole genome

# load libraries :  these will need to have already been installed (plus their dependencies.
library(GenomicLayers)
library(BSgenome.Scerevisiae.UCSC.sacCer3)    # using a "small" genome
genome <- BSgenome.Scerevisiae.UCSC.sacCer3

# set up a LayerSet on the genome :  a list with link to genome and GRanges objects to store position information
layerFullGenome <- createLayerSet.BSgenome(genome=genome, 
                                         layer.names=c("CpG_island",  "H3K27me1",  "H3K27me2","H3K27me3"),
                                           n.layers=4)

# create one or more binding factors:  an entity (e.g. a transcription factor) that recognises and, optionally, marks a specifc part of the genome
CGI<- createBindingFactor.DNA_regexp("CGI", patternString="(CG.{0,4}){3}CG", patternLength=20,
                                     mod.layers = "CpG_island", mod.marks=1, stateWidth=20)

bf.K27me3 <- createBindingFactor.layer_region("bf.K27me3", type="layer_region", 
                                          patternLength = 1, patternString = "N",
                                          profile.layers = "CpG_island",
                                          stateWidth=150, profile.marks = 1, 
                                          mod.layers = "H3K27me3", mod.marks = 1)

# add all the binding factors to a list
bfSet <- list(CGI=CGI, bf.K27me3=bf.K27me3)   # binding factors must be in a list and named. Easiest to use each BFs name.


system.time(
  newLayerGenome <- runLayerBinding.BSgenome(layerList=layerFullGenome, factorSet=bfSet, iterations=100000) 
)
newLayerGenome$layerSet
