

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


# Load list of factors from supplemental file
hp.motifs <- read.delim("data/XieKellis2005Nature/XieKellis2005Nature_S2.tab")
head(hp.motifs)


# create unique and informative names. (can use column No. to make unique)
hp.motifs$factor_name <- paste("XieKellis2005_", hp.motifs$No. ,sep="")
hp.motifs$factor_name <- ifelse(hp.motifs$Known_factor == "-", hp.motifs$factor_name , paste(hp.motifs$factor_name, hp.motifs$Known_factor, sep="_"))

# convert 'n' to 'N'
hp.motifs$Discovered_motif <- gsub("n", "N", hp.motifs$Discovered_motif)


# want to create a factorSet 
# each binding factor like this:-
#  bindingFactor <- list(name=name, type=type, 
#                        profile=profileList, 
#                       mods=modList)

#motifLibrary <- list()

# BUT the matching profile and modList for each factor is not specified and will need to be generated 
#	depending on the layerSet/layerList in use.
# We want the names, the motifs and the position_bias.

hp.motifs$hitsInTarget <- 0
i <-1
for(i in 1:nrow(hp.motifs)) {
thisPattern <- hp.motifs$Discovered_motif[i]

	hits <- matchPattern(thisPattern , targetSeq , fixed=FALSE)
	print(paste(i,length(hits)))
	hp.motifs$hitsInTarget[i] <- length(hits)

}

hist(hp.motifs$hitsInTarget, breaks=50)

# compare with a random set (should be very bad)

random.patterns <- data.frame()
for(i in 1:nrow(hp.motifs)) {
	patternLength <- max(4,rnbinom(1, 50, mu= 12))  # patterns must be at least of length 1
      pattern <- paste(sample(names(IUPAC_CODE_MAP), patternLength, replace=T), collapse="")
	hits <- matchPattern(pattern , targetSeq , fixed=FALSE)

	thisRow <- data.frame(pattern=pattern, hits=length(hits))
	random.patterns <- rbind(random.patterns, thisRow)
}

?nchar
mean(nchar(as.character(random.patterns$pattern)))
mean(nchar(as.character(hp.motifs$Discovered_motif)))
# the random patterns are longer (mean 12.24 vs 9.9) but often have 1000s of hits in the test_sequence
# all the real-life (but 'discovered') motifs have <= 130 hits and most have fewer than 20
# Information content? The real motifs are mostly more precise, with few Ns, and redundant codes.
# Could specify the letter freqeencies to match? 
paste(hp.motifs$Discovered_motif, collapse="")

?letterFrequency
# the letterFrequencies are far from even. Could be useful in generating 'random' binding factors
IUPAC.prob <- letterFrequency(DNAString(paste(hp.motifs$Discovered_motif, collapse="")), letters=names(IUPAC_CODE_MAP), as.prob=TRUE)
# note absence of not-X types (V, H, D, B)


# also the random patterns are sometimes shorter
range(nchar(as.character(random.patterns$pattern)))
range(nchar(as.character(hp.motifs$Discovered_motif)))

# remake random.patterns with freqs from the 'real' motifs and at least 6 characters
# used trial and error on rnbinom to get motif dist with same range
random.patterns.2 <- data.frame()
for(i in 1:nrow(hp.motifs)) {
	patternLength <- max(6,rnbinom(1, 50, mu= 7))  # patterns must be at least of length 1
      pattern <- paste(sample(names(IUPAC_CODE_MAP), patternLength, prob=IUPAC.prob,replace=T), collapse="")
	hits <- matchPattern(pattern , targetSeq , fixed=FALSE)

	thisRow <- data.frame(pattern=pattern, hits=length(hits))
	random.patterns.2 <- rbind(random.patterns.2, thisRow)
}

range(nchar(as.character(random.patterns.2$pattern)))
range(nchar(as.character(hp.motifs$Discovered_motif)))
mean(nchar(as.character(hp.motifs$Discovered_motif)))
mean(nchar(as.character(random.patterns.2$pattern)))
mean(hp.motifs$hitsInTarget)
mean(random.patterns$hits)
mean(random.patterns.2$hits)  # reduced but still nowhere near real.
sum(hp.motifs$hitsInTarget == 0)  # very few real motifs have zero hits in the target  # 8
sum(random.patterns$hits == 0)  # more of the random set have zero hits.  # 24
sum(random.patterns.2$hits == 0) #16
# So, this selection of 'real' motifs are far more specific than my letter-frequency and length matches 'random' sample.



