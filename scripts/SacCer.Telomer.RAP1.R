library(GenomicLayers)

## TESTING
#library(devtools)
#setwd(
## END OF TESTING

# Sir2p and Sas2p opposingly regulate acetylation of yeast histone H4 lysine16 and spreading of heterochromatin
# https://www.nature.com/articles/ng1017z
# Abstract:  The Sir3 protein helps form telomeric heterochromatin by interacting with hypoacetylated histone H4 lysine 16 (H4–Lys16). The molecular nature of the heterochromatin boundary is still unknown. Here we show that the MYST-like acetyltransferase Sas2p is required for the acetylation (Ac) of H4–Lys16 in euchromatin. In a sas2Δ strain or a phenocopy Lys16Arg mutant, Sir3p spreads from roughly 3 kb to roughly 15 kb, causing hypoacetylation and repression of adjacent chromatin. We also found that disruption of Sir3p binding in a deacetylase-deficient Sir 2Δ strain can be suppressed by sas2Δ. These data indicate that opposing effects of Sir2p and Sas2p on acetylation of H4–Lys16 maintain the boundary at telomeric heterochromatin.
# extract "By binding to cis-acting DNA sequences, Rap1p is believed to determine where heterochromatin initiates3. Silencing information regulators Sir2p, Sir3p and Sir4p are then recruited to Rap1p to form a complex that initiates and spreads along chromatin, interacting with the amino termini of histones H4 and H3 (ref. 4). These tails contain sites of acetylation on H4 at Lys5, 8, 12 and 16, and on H3 at Lys9, 14, 18, 23 and 27 (ref. 2). The interaction between Sir3p and Sir4p and the histone N termini probably involves hypoacetylated H4–Lys16 (but not H4–Lys5, 8 and 12), as determined by conservative and nonconservative amino-acid substitutions at these sites"

# Sir2p de-acetylates
# Sas2p acetylates H4K16ac

# https://www.science.org/doi/10.1126/science.2237406
# "Because RAP1 binds in vitro to the poly(C1-3A) repeats of telomeres, it has been suggested that RAP1 may be involved in telomere function in vivo. "

#https://www.sciencedirect.com/science/article/pii/S0955067497800117?via%3Dihub
# Grunstein (1997)  Review of molecular model
# Sir3 - increased dosage causes extension of silenced telomere ends to cover a gene.
# Suggest additional acetylation elsewhere in genome at potential Sir3 sites prevents 
#   Sir3 recruitment and therefore concentrates Sir3 at the telomeres.
#   See also note below about similar potential for H3k4me to have this effect

# See also Kimura et al (2002) with similar finding for Sas2p and Sir2p
# https://www.nature.com/articles/ng993z

# H3K4me at promoters may help keep Sir3p at telomeres as Sir3p cannot bind when H3K4me present
# https://www.jbc.org/article/S0021-9258(20)67880-2/fulltext
# Methylation of H3 Lysine 4 at Euchromatin Promotes Sir3p Association with Heterochromatin
# Santos-Rosa et al. (2004) J. Biol. Chem.

# Binding factors  ===================
# RAP1  (RAP1P?)  to bind to regexp
# patternLength set to minimum length the pattern can match (CA)
bf.RAP1 <- createBindingFactor.DNA_regexp("bf.RAP1", patternString="(C{1,3})A", patternLength=2,
                                     mod.layers = "Sir3p_potential", mod.marks=1, stateWidth=20)



# Sir2p - an HDAC. Is it specific to H4K16ac?
bf.Sir2p <- createBindingFactor.layer_region("bf.Sir2p", type="layer_region", 
                                              patternLength = 10, patternString = "N",
                                              profile.layers = "H4K16ac",
                                              stateWidth=10, profile.marks = 1, 
                                              mod.layers = "H4K16ac", mod.marks = 0)

# Sas2p acetylates H4K16ac
bf.Sas2p <- createBindingFactor.layer_region("bf.Sas2p", type="layer_region", 
                                             patternLength = 10, patternString = "N",
                                             profile.layers = "H4K16ac",
                                             stateWidth=10, profile.marks = 0, 
                                             mod.layers = "H4K16ac", mod.marks = 1)

# Sir3p
# binds to telomeres and is also spread 
# cannot bind acetylated regions or H3K4me
bf.Sir3p <- createBindingFactor.layer_region("bf.Sir3p", type="layer_region", 
                                             patternLength = 10, patternString = "N",
                                             profile.layers = c("H3K4me", "H4K16ac", "Sir3p_potential"),
                                             stateWidth=10, profile.marks = c(0,0,1), 
                                             mod.layers = "bound.Sir3p", mod.marks = 1)
# need also something to spread Sir3p along the chromosome.
# want a binding factor that creates mods either up or downstream of binding site. 
# New parameter 'offset.method' to accept function as argument
upDownFunc <- function(x)  {
  return(sample(c(x, -x), 1))
  
}
# PROBLEM cannot spread into region that has H3K4me or H4K16ac
bf.Sir3p.spread <-  createBindingFactor.layer_region("bf.Sir3p.spread", type="layer_region", 
                                            patternLength = 10, patternString = "N",
                                            profile.layers = c("H3K4me", "H4K16ac", "bound.Sir3p"),
                                            stateWidth=10, profile.marks = c(0,0,1), 
                                            mod.layers = "bound.Sir3p", mod.marks = 1, offset = 15,
                                            offset.method=upDownFunc)
# need to check these parameters to make sure it can't activate across an insulator

# this is just a hunch that I need something to displace Sir3p.
bf.eject.Sir3p <- createBindingFactor.layer_region("bf.eject.Sir3p", type="layer_region", 
                                                   patternLength = 10, patternString = "N",
                                                   profile.layers = "bound.Sir3p",
                                                   stateWidth=10, profile.marks = 1, 
                                                   mod.layers = "bound.Sir3p", mod.marks = 0)


# combine all bfs into a list to pass to runLayerBinding. Order DOES matter.


bfSet <- list(bf.RAP1=bf.RAP1, 
              bf.Sir2p=bf.Sir2p,
              bf.Sas2p=bf.Sas2p,
              bf.Sir3p=bf.Sir3p,
              bf.Sir3p.spread=bf.Sir3p.spread,
              bf.eject.Sir3p=bf.eject.Sir3p
              )   # binding factors must be in a list and named. Easiest to use each BFs name.
# specify abundances of bfs.
tf.abs <- structure(c(200, 200, 200, 200, 1000, 100), names=c(names(bfSet)))



# Layers needed: -
#  H4K16ac (possibly other aceytylations)  - hypo-acetylation required at telomeres
# H3K4me  - stops binding of Sirp3
# RAP1 binding/recruitment - can act as recruiter for other repressive modifiers.
# bound.Sir3p
library(BSgenome.Scerevisiae.UCSC.sacCer3) 

genome <- BSgenome.Scerevisiae.UCSC.sacCer3

# hack from https://support.bioconductor.org/p/83588/  to keep only some chromosomes within a BSgenome object.
keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}


sequences_to_keep <- "chrVII"   # telomere has many C{1,3}A repeats
genomeSub <- keepBSgenomeSequences(genome, sequences_to_keep)
genomeSub    # this should now still be a useable BSgenome object but with only one chromosome.  

# set up a LayerSet on the genome :  a list with link to genome and GRanges objects to store position information
layerSubGenome <- createLayerSet.BSgenome(genome=genomeSub, 
                                           layer.names=c("Sir3p_potential",  "bound.Sir3p",  "H3K4me","H4K16ac"),
                                          n.layers=4)

# probably worth testing the RAP1 bf
matchBindingFactor.BSgenome(layerSet=layerSubGenome, bindingFactor = bf.RAP1)
matchBindingFactor.BSgenome(layerSet=layerSubGenome, bindingFactor = bf.Sas2p)

# Run the simulation once to test it works.
system.time(
  newLayerGenome <- runLayerBinding.BSgenome(layerList=layerSubGenome, factorSet=bfSet, bf.abundances = tf.abs) 
)


newLayerGenome
tf.abs <- structure(c(1000, 200, 200, 200, 400, 100), names=c(names(bfSet)))
max.Sirp <- tf.abs['bf.Sir3p']
# Run the simulation, decreasing the availability of Sirp3 each iteration (could be function of binding).
newLayerGenome <- layerSubGenome
for(i in 1:500) {
  print(i)
  newLayerGenome <- runLayerBinding.BSgenome(layerList=newLayerGenome, factorSet=bfSet, bf.abundances = tf.abs, verbose=TRUE, collect.stats=TRUE) 
  # quick way of generating simple "movie". Only works when sim is one chrom.
  hist(start(newLayerGenome$layerSet$bound.Sir3p), breaks=20000, xlim=c(0, 100000), main=i, xlab="chr position")
  tf.abs['bf.Sir3p'] <- max(0, (max.Sirp - length(newLayerGenome$layerSet$bound.Sir3p)) )  # reduce availability of bound Sir3p
}

hist(start(newLayerGenome$layerSet$bound.Sir3p), breaks=20000,  main=i, xlab="chr position")

hist(start(newLayerGenome$layerSet$bound.Sir3p), breaks=1000,  main=i, xlab="chr position")

plot(newLayerGenome$history[,"Coverage.Sir3p_potential"])
plot(newLayerGenome$history[,"Coverage.bound.Sir3p"])   # interesting. Semi-linear.
# ideally, need chrom distribution.

plot(x=start(newLayerGenome$layerSet$bound.Sir3p),y=rep(1, length(newLayerGenome$layerSet$bound.Sir3p)), xlim=c(1,seqlengths(genomeSub)["chrVII"]), xlab="Chrom position", ylab="Sir3p bound")


rap1.hits <- matchBindingFactor.BSgenome(layerSet=layerSubGenome, bindingFactor = bf.RAP1)
hist(start(rap1.hits), breaks=100)  # seems fairly uniform across chromosome
hist(start(rap1.hits), breaks=20000, xlim=c(0, 1000))



#hist(width(newLayerGenome$cache$bf.RAP1$LAYER.0), breaks =40)  # the 30+ hits at the telomere are serious outliers and should act lile magnets for RAP1
hist(start(newLayerGenome$layerSet$bound.Sir3p), breaks=100)
hist(start(newLayerGenome$layerSet$bound.Sir3p), breaks=20000, xlim=c(0, 1000))

# test to see what happens with just the Sas2p - should add lots of acetylation

bfSet.Ac <- list(bf.Sas2p=bf.Sas2p)
newLayerGenomeAc <- runLayerBinding.BSgenome(layerList=layerSubGenome, factorSet=bfSet.Ac, bf.abundances = 1000, verbose=TRUE, collect.stats=TRUE) 

genomeSub[[newLayerGenome$layerSet$bound.Sir3p[1]]]
getSeq(genomeSub, newLayerGenome$layerSet$bound.Sir3p)
hist(width(newLayerGenome$layerSet$bound.Sir3p), breaks=100)  # lots of frags that need to be removed?

## OBSERVATIONS -------------------------------

# Simulations with unlimited bf abundances cover the whole chromosome.
# Adding feedback on quantity of Sir3p so that available decreases as more is bound, 
# leads to punctate spikes of binding at 10-12kb distance
#    (weirdly evenly spaced: 0, 13, 20, 30, 40, 50, 60, 80, 100).
#     A second run produced similar spacing but at different locations so driven by random finding of initiator.
# Driven in part by high ratio of spreader (bf.Sir3p.spread) to initiator (bf.RAP1)
# > tf.abs
# bf.RAP1        bf.Sir2p        bf.Sas2p        bf.Sir3p bf.Sir3p.spread  bf.eject.Sir3p 
# 200             200             200               200*            1000             100 
#  N.B. possible ERROR the blocks/coverage suggest there are more than 200 bf.Sas2p bound.
# This is partly a function of the cleaner function and perhaps need to widen the kick off pattern
#   or run a cleaner so that items less than accepted width should be removed.

###DEBUGGING ---------------------------------
#First results with no H3K4me had no apparent focus to telomeres.
# Could be due to lack of H3K4me or could be because of way overlapping/consecutive DNA based hits
#    are reduced to one range and therefore not selected in response to frequency.

thisPattern <- "(C{1,3})A"
bsParams <- new("BSParams", X=genomeSub, FUN=gregexpr)  # set up params for using bsapply
grepResultBS <- bsapply(bsParams, pattern=thisPattern)  # run gregexpr over each chromosome separately
chromName <- "chrVII"
grepResult <- grepResultBS[[chromName]]
win.hits <- GRanges(chromName, IRanges(start= as.integer(grepResult[[1]]), width = attr(grepResult[[1]], which="match.length", exact=TRUE)),seqinfo=seqinfo(genome))
reduce(win.hits, ignore.strand=TRUE)   # this is causing loss of individual hit locations for regular expression. Seems a bad idea for regexp, or any DNA based hits.

# probably worth testing the RAP1 bf
matchBindingFactor.BSgenome(layerSet=layerSubGenome, bindingFactor = bf.RAP1)  # current implementation reduces adjacent hits to one, which is bad.
GenomicLayers::matchBindingFactor.BSgenome(layerSet=layerSubGenome, bindingFactor = bf.RAP1) 

# If Sir3p stays bound to telomeres, then it's abundance may decline as simulation runs. 
# Ultimately, this would reduce it's spread from the telomere and the turnover may 
# explain the eventual pattern.  
# Probably needs a Low-frequeny way to remove this from the chromosome. 
# Grun