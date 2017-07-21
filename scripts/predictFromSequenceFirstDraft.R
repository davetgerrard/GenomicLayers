

# test script to predict events from sequence with only a sequence of motif-like binding events.


# data should be genome sequence, or set of long sequences
# On top of this, n layers (1-100?) (start with 5). 
# That, along with the sequence, can be read by the binding factors
# Factors bind (probabilisticall) where they match their motif and then alter specific layers
# (the sequence is immutable)




# LAYERS:  matrix or list of vectors as long as the DNA sequence. One 'row' per 'layer'
# Represent nucleosome binding, histone mods, methylation, compaction, TF binding etc.
# 


# FACTORs: anything that recognises the DNA + LAYERS, binds to it and alters the state of the layers
#  each FACTOR has a multi-dimensional motif 
#  each FACTOR has an EFFECT
#  each EFFECT has a target LAYER, an AREA_OF_EFFECT and a EFFECT_DISTANCE. 
#  The area and distance may be sampled from a distribution

# e.g. FACTOR_X  binds the pattern ACACCAAC IF LAYER-2 is in state 1 
# or more sophisticated a MATCH could be the intersection of matches on the set of layers
# i.e. work out matches on sequence, work out matches on all layers (where relevant), 
# then full factor matches are the intersection of all the layer matches. 
#  LAYER-0 (SEQUENCE):   TCGATATATACGCCGCTATACGGCTTAGCT
#  LAYER-1 :             111111111111110000000000000000
#
#  MATCH-0 (TATA):       000011111100000011110000000000     # N.B. two overlapping matches....
#  MATCH-1 (1)           111111111111110000000000000000
#
# FULL FACTOR MATCH:     000011111100000000000000000000
#                            ^^^^^^



# ?N.B.  The overlapping matches above suggest it might be better to record each distinct match in a list of match positions
# e.g. MATCH-0 (TATA):  5-8,7-10,17-20


# MOTIF-structure
# FACTOR
#   |----> MOTIF
#           |----> LAYER-0
#           |----> LAYER-1
#           |----> LAYER-n



# FACTORS should be available at different frequencies (different abundancies).
# In the first case, the abundancies can be constant. (house-keeping, steady state).
# Later can be developed to alter during 'development'

# FACTORS should also be applied with a (nearly) determined order.
# I'm not sure how important this is. But am thinking it will mirror the biological state of the sequence of binding events.
# Perhaps for house-keeping genes, a system could be evolved for which the order is not important and the result is just
# an outcome of the availability (abundance) of factors and the sequence.  
# This isitself an interesting biological question: how much is the order of regulation important?
#  (or more specifically, 

# How many ROUNDS of application should occur? 

# What is the initial state of the DNA + LAYERS. Sequence plus all zero?



# How to 'evolve' the factors.
# Motif for each layer (mostly 1s and 0s) but special for DNA-motif
# pattern, 
# Effect of each factor:  layer(s), range/distance/area
# e.g. 
# bind at position 0, set  LAYER-3 to 1 @ +5-+20
# bind at position 0, set LAYER-
# N.B. LAYER-0 (DNA) cannot be edited.
# Design:  for each layer/motif within each factor, it might be good to give it a MUTABLE flag.
# i.e. we might want to force some factors to only operate on certain layers (e.g. the final assesment factor could be ties to an assessment layer)


# Half a thought.  The search space for motifs is very large (see note below on nucleosomes).
# Could FACTORS be down-weighted if they never even match to the DNA?


# The final assessment will be a ability of one of the factors (any or a pre-defined one?) to match
# the multi-layer profile in a set of positions 'similar' to a training set.
# e.g. overlap a set of house-keeping TSSs.
# How to measure agreement?   (overlap within a window?)





# Future developments: 
#   seed with motifs from JASPAR, Jolma et al., etc.
#   tissue- or stage-specific genes  - need to link to changing abundance of factors.








# Problems/inadequacies.
# Fails to model interactions
# Fails to model sliding objects (e.g. modified histones within nucleosomes)
#     nucleosomes could be modelled if they recognise and bind to regions with gap larger than themselves. 
#     e.g. motif is gap (state 0) of 200bp, the 'factor' binds in the centre and sets the same layer to filled (1) with size ~140bp





# Coding notes


#http://stackoverflow.com/questions/14249562/find-location-of-character-in-string
library(stringr)
str_locate_all(pattern ='2',"the2quickbrownfoxeswere2tired")
str_locate_all(pattern ='TATA',"TCGATATATACGCCGCTATACGGCTTAGCT")  # does not list all overlapping patterns.
str_locate_all(pattern ='TATA',"TATATATATATATATATATA") 
gregexpr(pattern ='2',"the2quickbrownfoxeswere2tired")
gregexpr(pattern ='TATA',"TCGATATATACGCCGCTATACGGCTTAGCT")
gregexpr(pattern ='TATA',"TATATATATATATATATATA")
# How to find ALL matches even over-lapping ones
#http://stackoverflow.com/questions/14863026/javascript-regex-find-all-possible-matches-even-in-already-captured-matches/14863268#14863268

# How about dedicated motif finding R packages?
#http://www.niravmalani.org/finding-motif-of-interest-in-a-genome-using-r/
library(Biostrings)
hits1 <- matchPattern(DNAString("TATA"), DNAString("CGATATATACGCCGCTATACGGCTTAGCT"), fixed = FALSE)  # DOES FIND ALL OVERLAPPING MATCHES.
hitsN <- matchPattern(DNAString("N"), DNAString("CGATATATACGCCGCTATACGGCTTAGCT"), fixed = FALSE) 
# note that BString object may be useful.
hits2 <- matchPattern(BString("1111"),BString("000000111110000001111"))  # but not sure how to do regular expressions with this
# e.g. want to match "0(1-100)1(1-100)0(1-100)"
matchPattern(BString("1111"),BString("110100000111110000001111"), max.mismatch=1)
matchPattern(BString("1111"),BString("110100000111110000001111"), max.mismatch=1,with.indels=T)   # only the best local match.
# How to do really long matches?  
# Might be simplist just to have long string with approximate pattern but allow quite a few mismatches.
#e.g. 
chrom <- BString("000000000000000000000110000000000000000000000111111111111111111111111111111111111111100000000000000000000000000000000000000000000000000000000001000000000000000000001111111111111111111011111111111111111111110000000000000000000000000000000000000000000110000000000000000000001000000000000011111111111111111111111111110111111111111100000000000000000000000000000000000000000000000000000000000000000000000000000001111111111111111111111111111111111111111110000000000000000000000000000000000000000000")
pattern <- BString("000000000000000000000000000111111111111111111111111111111111111111111000000000000000000000000000000000000000")
mismatch.rate <- .1
matchPattern(pattern,chrom , max.mismatch=round(length(pattern)*mismatch.rate))
# one slight problem with this might be that each block match produces very many similar approximate matches. 
#     But, if it is a case of a fACTOr binding the dna and doing something, it may not matter 

# see also matchPDict



# layer marking test
# Worried that the biostrings objects may not be easy to modify

require(Biostrings)
#targetSeq <- DNAString('TGCGTTGAGC')
targetSeq <- readDNAStringSet("C:/Users/Dave/Dropbox/Temp/hg19.HOXA.fa")[[1]]
emptyLayer <-  BString(paste(rep(0,length(targetSeq)), collapse=""))

# make a slit of layers, with the target sequence as LAYER.0
layerSet <- list(LAYER.0 = targetSeq)
layerSet[['LAYER.1']] <- emptyLayer

# can a layer be edited?
test <- emptyLayer
test[1] <- "1"
# Yes


# define a binding factor and it's influence
bindingFactor <- list(name="TestFactor", 
                      profile=list(LAYER.0=list(pattern=DNAString("ATCCNCCTTCCTT") , mismatch.rate=0), 
                                   LAYER.1=list(pattern=BString("000000000000000"), mismatch.rate=0)), 
                      mods=list(LAYER.1=list(state="1", stateWidth=20, offset=0, align="centre")))







# might be a saving by not calculating every hit when a layer is completely uniform or not tested.
# currently DNA matches using fixed length string with IUPAC codes (adapt later for pwm matching)

theseHits <- matchBindingFactor(layerSet, bindingFactor)




modifyLayer(layerSet, 10:20, layer="LAYER.1")
# could this be written as function on an object?
# e.g. 

# for each valid hit, modify layerSet (allowing for multiple mods across layers for each hit)
# for each hit, need a single coordinate for the modifications to be based upon
# for now, use midpoint
thisHit <- theseHits[1]
thisHitPosition <- start(thisHit) + floor(width(thisHit)/2)


modifyLayerByBindingFactor(layerSet, position=thisHitPosition, bindingFactor=bindingFactor)



# next is to create more bindingFactors, and perhaps more layers and then iteratively perform many layerMods. 
# Also need to build imprecise matching and tolerances for such. will require additional level under 'profile'
#     so that each layer match can have a different tolerance.
# Some example binding Factors:-
#   imperfect matching of many GCs?
IUPAC_CODE_MAP


allFactors <- list()

allFactors[['GC.rich']] <- list(name="GC.rich", 
                                profile=list(LAYER.0=list(pattern=DNAString("SSSSSSSSSSSSSSSSSSSSSSS"), mismatch.rate=0.1) ,
                                             LAYER.1=list(pattern=BString("000000000000000"), mismatch.rate=0)), 
                                mods=list(LAYER.1=list(state="1", stateWidth=20, offset=0, align="centre")))

allFactors[['TATA']] <- list(name="TATA", 
                                profile=list(LAYER.0=list(pattern=DNAString("TATAA"), mismatch.rate=0) ,
                                             LAYER.1=list(pattern=BString("000000000000000"), mismatch.rate=0)), 
                                mods=list(LAYER.1=list(state="1", stateWidth=10, offset=0, align="centre")))

# cleaner bindingFactor removes marks from layer.1
allFactors[['cleaner.LAYER.1.width.1']] <- list(name="cleaner.LAYER.1.width.1", 
                             profile=list(LAYER.0=list(pattern=DNAString("N"), mismatch.rate=0) ,
                                          LAYER.1=list(pattern=BString("1"), mismatch.rate=0)), 
                             mods=list(LAYER.1=list(state="0", stateWidth=1, offset=0, align="centre")))


theseHits <- matchBindingFactor(layerSet, allFactors[['GC.rich']])

theseHits <- matchBindingFactor(layerSet, allFactors[['TATA']])
length(theseHits)
theseHits <- matchBindingFactor(layerSet, allFactors[['cleaner.LAYER.1.width.1']])   # no hits on clean layer
length(theseHits)
# first modify layer.1, then search again.
moddedLayerSet <- modifyLayerByBindingFactor(layerSet, position=thisHitPosition, bindingFactor=bindingFactor)
theseHits <- matchBindingFactor(moddedLayerSet, allFactors[['cleaner.LAYER.1.width.1']])
modifyLayerByBindingFactor(moddedLayerSet, position=5, bindingFactor=allFactors[['cleaner.LAYER.1.width.1']])  ## removes mark on layer.1 at given position.



bindingFactor <- list(name="TestFactor", 
                      profile=list(LAYER.0=DNAString("ATCCNCCTTCCTT") , LAYER.1=BString("000000000000000")), 
                      mods=list(LAYER.1=list(state="1", stateWidth=20, offset=0, align="centre")))




## need to run multiple iterations of binding against the same sequence.

bindFreqs <- rep(1, length(allFactors))           # used to store the (relative) numbers of binding factors, here all equal
sample(names(allFactors), size=30,prob=bindFreqs, replace=T)   # generate random vector of binding factor names







modLayerSet <- runLayerBinding(layerSet=layerSet, factorSet = allFactors, iterations=100)

letterFrequency(modLayerSet$LAYER.1, letters=c("0", "1"))   # how many of layer.1 were set to 1

modLayerSet <- runLayerBinding(layerSet=layerSet, factorSet = allFactors, iterations=1000)






# could do with some plots to show where on a sequence layers have been changed
#as.matrix(modLayerSet$LAYER.1[1:100])
layer.vector <- as.numeric(strsplit(as.character(modLayerSet$LAYER.1),"")[[1]])
plot(layer.vector, type="l")  # very slow.

# quick test to see if the same 'barcode' is created after two runs 
par(mfrow=c(2,1))
modLayerSet1 <- runLayerBinding(layerSet=layerSet, factorSet = allFactors, iterations=1000)
layer.vector.1 <- as.numeric(strsplit(as.character(modLayerSet1$LAYER.1),"")[[1]])
plot(layer.vector.1, type="l")  # very slow.
modLayerSet2 <- runLayerBinding(layerSet=layerSet, factorSet = allFactors, iterations=1000)
layer.vector.2 <- as.numeric(strsplit(as.character(modLayerSet2$LAYER.1),"")[[1]])
plot(layer.vector.2, type="l")  # very slow.


cor(layer.vector.1, layer.vector.2)   # first try gave correlation of 0.84 between two runs with 1000 iterations each.
# suggests that the majority of bindable sites have been found 
# Are they at a steady state?




# Need to test resulting pattern in a target layer against known set of regions
# e.g. positions of transcriptions start sites (TSS)
# for this need real data set with real TSS locations.
transcriptTable <- read.delim("C:/Users/Dave/Dropbox/Temp/hg19.HOXA.transcripts.tab")
head(transcriptTable)
base.0 <- 27106375   # one before the start co-ordinate used to get this DNA sequence.

tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.

tss.vector <- rep(0, nchar(targetSeq))
tss.vector[tss.positions] <- 1
plot(tss.vector, xlim=c(0, nchar(targetSeq)), type="l")
# need a function(s) to generate bindingFactors
# plot tss locations, then the layer.1 from a few patterns
par(mfrow=c(3,1))
plot(tss.vector, xlim=c(0, nchar(targetSeq)), type="l")
plot(layer.vector.1, type="l")
plot(layer.vector.2, type="l")

cor(tss.vector, layer.vector.1)   # no correlation, not a problem for the first run - something to work on !
cor(tss.vector, layer.vector.2)  




  
  
#  IUPAC_CODE_MAP
  


createRandomBindingFactor("testRandom", layerSet, type="DNA_motif", test.layer0.binding=FALSE, test.mismatch.rate=.1 ) 


pattern.length <- 20
paste(sample(names(IUPAC_CODE_MAP), size=pattern.length, replace=T), collapse="")  # uniform probabilities.
# how to incorporate homogeneity? or is that something to evolve?

# Strategies for 'evolving' improved TSS predictors
# 1. seed a single factorSet, then mutate it for improvements
# 2. seed many factorsets, test them all and find the best. re-run it for stability.

# Going to try option 2 to start with, just to test mechanics.
# need many factorSets.
# Need a layerSet with more than 2 layers.
# need to specify a target layer and assessment criteria. (in this case can do that at the end)

n.layers <- 5

layerSet.5 <- list(LAYER.0 = targetSeq)
for(i in 1:n.layers) {
  layerSet.5[[paste('LAYER.', i , sep="")]] <- emptyLayer
}

factorSetRandom <- list()
for(i in 1:50) {
  factorSetRandom[[paste("bf.",i ,sep="")]] <- createRandomBindingFactor(paste("bf.",i ,sep=""), layerSet.5, type="DNA_motif", test.layer0.binding=FALSE, test.mismatch.rate=.1 ) 
  
}

print.bfSet(factorSetRandom)


modLayerSet <- runLayerBinding(layerSet=layerSet.5, factorSet = factorSetRandom, iterations=100)
lapply(modLayerSet[-1], letterFrequency, letters=c("0", "1")) 


theseHits <- matchBindingFactor(layerSet.5,  factorSetRandom[[1]])
moddedLayerSet <- modifyLayerByBindingFactor(layerSet.5, position=thisHitPosition, bindingFactor= factorSetRandom[[1]])






# would be good to test novel binding factors with DNA motifs to see if they ever exist in
# target sequences?





# later on, when it comes to mutations, 
# distinct functions for mutating integers, fractions, and binding profiles (DNA and binary)
# have separate occurrence frequencies for each type of mutation. 
# chance of mutation of an individual binding factor should just be 1/n.factors * mutation.rate
# so 1. calc number of new mutations
#     2. calc type of mutations (can be based on weighted probs)
#     3. pick a bindingFactor at random
#     4. mutat appropriately, allowing 1 <-> 1 muts (e.g. unchanged).
#  Is it fair to compare mutations in large region binding with small motifs? - does it matter?
# Might be able to work around homogeneity issue but having special class of homgeneous indels
#   where a letter is selected from the existing pattern and an insert made next to if of the same type
#   But then, how to do homogeneous deletion - only valid for patterns that have homogeneous base clusters.
#      Could be measured by RunLengthEncoding (rle) ? 


# Need to factor in strand. Perhaps by search for all patterns on both strands?
# only layer.0 (the DNA-seq) will be reverse-complemented
# other layer matches will just be reversed


