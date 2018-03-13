
# Pluripotency model 1 ------------------

# Pre-pattern the genome with NANOG, OCT4, SOX2 etc.

# Then allow establishment of transcription at these sites.

# Also allow slow increase in polycomb at non-marked sites.

# Count genes as active or inactive at stages through the simulation and compare with patterns of polycomb repression in early embryos.

library(BSgenome.Mmusculus.UCSC.mm9)
library(GenomicLayers)


# some motif files (from HOMER) in ../data/homerMotifs/motifs/

# e,.g. Oct4 'oct4.motif'

# >ATTTGCATAA	Oct4(POU,Homeobox)/mES-Oct4-ChIP-Seq(GSE11431)/Homer	7.325577	-5.831779e+03	0	55120.0,27560.0,7917.0,6798.0,0.00e+00
# 0.852	0.048	0.048	0.052
# 0.001	0.001	0.001	0.997
# 0.001	0.043	0.001	0.955
# 0.155	0.001	0.001	0.843
# 0.107	0.001	0.891	0.001
# 0.008	0.990	0.001	0.001
# 0.997	0.001	0.001	0.001
# 0.001	0.004	0.001	0.994
# 0.676	0.047	0.125	0.152
# 0.434	0.107	0.063	0.396


# motifs would be good, but perhaps lets just go with consensus for now. 
# not yet worked out import from homer.motif format. see GenomicLayers/scripts/pwm.import.R



# define the genome to use 
genome <- BSgenome.Mmusculus.UCSC.mm9   # for convenience
genome
seqnames(genome)
summary(genome)
organism(genome)
provider(genome) 





thisChrom <- genome[["chrX"]] 


# set up a layerSet on chrom X.
print("set up layerSet")
layerSet.X <- list(LAYER.0 = thisChrom)
layerSet.X[['Active']] <- IRanges()
layerSet.X[['CpG_island']] <- IRanges()
layerSet.X[['PRC']] <- IRanges()
layerSet.X[['H3K27me3']] <- IRanges()

layerList.X <- list(layerSet=layerSet.X, history=NULL)  # add some metadata

bf.Oct4 <- createBindingFactor.DNA_motif("Oct4", patternString="ATTTGCATAA",
                                         profile.layers=NULL,profile.marks=NULL,
                                                mod.layers = "Active", mod.marks=1, stateWidth=200)

testBinding <- matchBindingFactor(layerSet = layerSet.X, bindingFactor = bf.Oct4,max.window = 10000000, verbose = TRUE)
length(testBinding)
result <- runLayerBinding(layerList=layerList.X, factorSet=list(bf.Oct4=bf.Oct4), iterations = 200)
# TODO allow passing of max.window to matchBindingFactor from within runLayerBinidng



# whole genome attempt...
mmLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE, layer.names=c(paste0("LAYER.", 1:4), "Active"))
#layerList.mm9 <- list(layerSet=mmLayerSet, history=NULL)
system.time(
  wgTest <- runLayerBinding.BSgenome(layerList=mmLayerSet, factorSet=list(bf.Oct4=bf.Oct4), iterations = 2000)
  , gcFirst = TRUE)

system.time(
  wgTest2 <- runLayerBinding.BSgenome(layerList=mmLayerSet, factorSet=list(bf.Oct4=bf.Oct4), max.window=100000000,iterations = 2000, , watch.function=function(x, ...){print("")})    # the watch.function required to get past problem in runLayerBinding.BSgenome()
  , gcFirst = TRUE)   # might be faster if can pass 'max.window' param to matchBindingFactor.BSgenome() - Nope, either not working or not helpful..

system.time(
  wgTest2 <- runLayerBinding.BSgenome(layerList=mmLayerSet, factorSet=list(bf.Oct4=bf.Oct4), verbose=T,max.window=100000000,iterations = 2000, , watch.function=function(x, ...){print("")})    # the watch.function required to get past problem in runLayerBinding.BSgenome()
  , gcFirst = TRUE)   # might be faster if can pass 'max.window' param to matchBindingFactor.BSgenome() - Nope, either not working or not helpful..


sum(width(wgTest$layerSet$Active))
sum(width(wgTest2$layerSet$Active))    # similar but not the same.
stopifnot(FALSE)
# 
# # genome-wide layerSets still VERY slow for DNA pattern matching.
# # use single chromosome to optimise and perhaps implement caching of genome-wide DNA matches.
# 
# # Create genome layerset
# source("../GenomicLayers/R/createLayerSet.BSgenome.R")  # not yet part of package?
# mmLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)
# 
# 
# #bf.CpGisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){9}CG", patternLength=20,
#  #                                              mod.layers = "CpG_island", mod.marks=1, stateWidth=200)
# 
# bf.Oct4 <- createBindingFactor.DNA_motif("Oct4", patternString="ATTTGCATAA",
#                                                mod.layers = "LAYER.1", mod.marks=1, stateWidth=200)
# 
# 
# result <- matchBindingFactor.BSgenome(mmLayerSet, bindingFactor = bf.Oct4, verbose = TRUE)  # VERY SLOW, AND WRONG! !
# 
# result <- runLayerBinding.BSgenome(layerList = mmLayerSet, factorSet=list(  bf.Oct4 = bf.Oct4),  verbose = TRUE)  # VERY SLOW, AND WRONG! !
# 
# 
# #?createLayerList.DNAstring
# #thisChrom <- genome[["chr1"]]
# #mmLayerSet <- createLayerSet.BSgenome(genome=thisChrom, n.layers = 5, verbose=TRUE)
# 

