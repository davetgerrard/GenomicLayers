library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

genome <- BSgenome.Hsapiens.UCSC.hg19
chrI <- genome$chrI
p1 <- "ACCCAGGGC"
countPattern(p1, chrI, max.mismatch=1)

seqnames(genome)
rcpattern <- reverseComplement(pattern)


# best to load one chromosome at a time, scan for all patterns, keep results and then log to next pattern.

thisChromName <- seqnames(genome)[1]
subject <- genome[[thisChromName ]]

p1 <- "ACCCAGGGC"
countPattern(p1, subject , max.mismatch=1)
#countPattern()


# how many possible patterns are there of a given length of DNA?
# using just the four pases

library(gtools)
combinations(3,2,letters[1:3])
combinations(3,2,letters[1:3],repeats=TRUE)

permutations(3,2,letters[1:3])
permutations(3,2,letters[1:3],repeats=TRUE)

#permutations


nrow(permutations(4,2,c("A","C","G","T"),repeats=TRUE))

for(i in 2:10)  {

print(paste(i,":", nrow(permutations(4,i,c("A","C","G","T"),repeats=TRUE))))
}

# there are too many possibilities!  To list them all.

subject <- genome[[thisChromName ]]
permuteTable <- permutations(4,2,c("A","C","G","T"),repeats=TRUE)
patterns <- character()
permuteCounts <- integer()
for(i in 1:nrow(permuteTable))  {
 patterns[i] <- paste(permuteTable[i,], collapse="")
 permuteCounts[i] <- countPattern(patterns[i], subject , max.mismatch=0)
}
permuteCounts
names(permuteCounts) <- patterns
barplot(permuteCounts)


# STRATEGIES  --------------

# 1.  Generated random patterns, determine their counts in the genome. 
# 2.  Parse the genome into patterns of length 'n'  
# 3. Use only known motifs.

# Option 1 is time consuming but will work with degenerate codes. May take forever to find long motifs.
# Option 2 is fast but not easy to then condense into degenerate patterns.
# Option 3 will include long degenerate motifs but will be constrained to well-studied factors, many with low number of binding events. Also will not find new motifs.



#  Homo-polymer patterns.
# for each member of IUPAC_code (except N?), determine number of hits of homo-polymer of length n
# differnet lengths and mismatches.

patternLength <- 8
pattern <- paste(rep(names(IUPAC_CODE_MAP)[8], patternLength), collapse="")

countPattern(pattern , subject , max.mismatch=0, fixed=FALSE)

mismatches <- 0
patternLength<- 8

patternCounts <- structure(integer( length=length(IUPAC_CODE_MAP)), names=names(IUPAC_CODE_MAP))
for(thisNuc in names(IUPAC_CODE_MAP))  {
	pattern <- paste(rep(thisNuc , patternLength), collapse="")
	print(pattern)
	patternCounts[thisNuc ] <- countPattern(pattern , subject , max.mismatch=mismatches, fixed=FALSE )
}

barplot(patternCounts, main=paste("length:", patternLength, "mismatches:", mismatches))




