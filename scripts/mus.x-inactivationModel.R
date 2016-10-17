

# Simulate 2-stage X-inactivation in mouse
#environment, packages and functions
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
library(BSgenome.Mmusculus.UCSC.mm10)
#require(Biostrings)   # included in mm10 package dependencies
source('scripts/pfs.functions.R')

# load X chromosome sequence 
genome <- BSgenome.Mmusculus.UCSC.mm10 
thisChrom <- genome[["chrX"]] 

# set up a layerSet on chrom X.

layerSet.X <- list(LAYER.0 = thisChrom)
layerSet.X[['CpG_island']] <- IRanges()
layerSet.X[['PRC']] <- IRanges()
layerSet.X[['H3K27me3']] <- IRanges()

layerList.X <- list(layerSet=layerSet.X, history=NULL)



# specify binding factors in model (and proportions?)

#   Stage 1.  - Gene rich regions  
#               CpG islands    [Pinter et al., 2012](http://www.citeulike.org/user/daveGerrard/article/11359142)
#   Stage 2. - Spread out from these regions in a two stage feedback (could start with one factor?).
#           - Meant to represent recruitment of PRC and then marking of H3K27me3.


bf.CpG <- createBindingFactor.DNA_motif("CpG", patternString="CG", profile.layers=NULL,profile.marks=NULL, mod.layers = "CpG_island", mod.marks=1)

bf.PRC <- createBindingFactor.layer_region("PRC",  patternLength=200, mismatch.rate=.4, profile.layers = "CpG_island", profile.marks = 1,  mod.layers = "PRC", mod.marks=1)

bf.CpGisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){30}CG", patternLength=300, profile.layers=NULL,profile.marks=NULL, mod.layers = "CpG_island", mod.marks=1)
#N.B. temporarily setting a fixed patternLength, because don't know how to implement variable patternLength yet.

# want a binding factor that creates mods either up or downstream of binding site. New parameter 'offset.method' to accept function as argument
upDownFunc <- function(x)  {
  return(sample(c(x, -x), 1))
  
}

bf.spreadRep <- createBindingFactor.layer_region("spreadRep", patternLength=200, mismatch.rate=.4, profile.layers = "PRC", profile.marks = 1,  mod.layers = "H3K27me3", mod.marks=1, stateWidth=500,offset=1000, offset.method=upDownFunc)



results1 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpG)    # should be many hits, each of 2bp.  
results2 <- matchBindingFactor(layerSet.X, bindingFactor = bf.PRC)    # should be no hits.
results3 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland)
results4 <- matchBindingFactor(layerSet.X, bindingFactor = bf.spreadRep)    # should be no hits.

###!!!!
# How to combine/detect clusters of high density?
# this kind of matching is not implemented!
#   the CpG factor is matching a random selection of 2bp CpG di-nucleotides.  To then do inexact matching of regions, is not obvious using iRanges intersections.

# could instead do it as a CpG island regexp?     /(CG.{0,9}){3}CG/    # finds at least 4 (3 plus 1) CGs with gap of 0 to 9 between each.
#     http://regexr.com/3edjk
# Raises the question of what the databases class as a CpG island (https://www.biostars.org/p/79046/) vS. what might be recognised by a binding factor.
#     Does a TF care about local over-representation?  

# Biostrings:::matchLRPatterns()
# Biostrings:::gregexpr2("aa", c("XaaaYaa", "a"))   # only works in a fixed mode, but would give locations of all CG quickly.

# implemented createBindingFactor.DNA_regexp()  and added to matchBindingFactor
#  BUT, sensing whole CpG islands is perhaps not how it works. Or is it?    
# could use CpGisland BF that must not have high % methylation in other layer
#     BUT, how to do overlaps with % matching using IRanges?


# runLayerbinding
#     Need to match both strands....


#   combine factors into factor set (list)
#XFS <- list(bf.CpG=bf.CpG, bf.PRC=bf.PRC, bf.CpGisland =bf.CpGisland ,bf.spreadRep=bf.spreadRep)
XFS <- list(bf.CpGisland =bf.CpGisland , bf.PRC=bf.PRC, bf.spreadRep=bf.spreadRep)

mod.X <- runLayerBinding(layerList.X, factorSet = XFS, iterations = 10000, verbose=T)

# TODO would be good to add turning off of CpG islands as repression spreads. Representing methylation. Could be modelled as separate layer or 


# compare with x-inactivation data.


#optimise binding factors to better match sequential inactivation



 biocLite("BSgenome.Mmusculus.UCSC.mm10")
 