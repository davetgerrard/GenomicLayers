
# Do CpG islands stand out using simple grep on a patten of N CpGs with up to G nucleotides between them?
# short answer: no.     
# off shoot from mus.x-inactivationModel.R
# Raises the question of what the databases class as a CpG island (https://www.biostars.org/p/79046/) vS. what might be recognised by a binding factor.
#     Does a TF care about local over-representation?  

setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
library(BSgenome.Mmusculus.UCSC.mm10)
#require(Biostrings)   # included in mm10 package dependencies

# load X chromosome sequence 
genome <- BSgenome.Mmusculus.UCSC.mm10 
thisChrom <- genome[["chrX"]] 


testString <- "CATATGCGCTGCTCTAATAGACTCGCCGCGCTAGCTACGCGTCACGTTAGCTATGCTACGTTGCTGACTAGCTCGATCGTGCTATGAATATATCGTCGATGCTGAGCGCTCTCCTCGCGCCGGCGGCGGGGCGCCCGCCGCGCCGCGTATATATATGCATGCCGCGTATATAGCGCGCCGTCCTGAATATATCGCGCTCTAATATATATCGGCGCTTAATTGATATAATCGCCGTGATCTAGCTAGTCATGTCGAGGACAGTCAGTAGAGATCCATGCTGTTGTGCTCTCTCTCTAGAGGAGAGACCAGAGTCTCTCTCCGGGCTGAAGAGTATGCATATGCATGCTAGCATGCATGCATGCATGCTAGCATGCATGCATGCTGCGCAGTCATGCATGACGTCAGTACCGTGCGACTGAGCTACGTATCAGATCCTGATCACGCGCGATTAATTGCTCTAGTCTGATATGCTGACGCGCGCGATTTACTTAATGGTTATATTGGGGATTCTTCCCCAATG"

gregexpr("(CG.{0,9}){3}CG", testString)
# but would this work on a chromosome
testResult <- gregexpr("(CG.{0,20}){30}CG", thisChrom)
# (yes it does!)
# how wide are they
hist(attr(testResult[[1]], which="match.length", exact=TRUE))   # how to access widths.
median(attr(testResult[[1]], which="match.length", exact=TRUE)) 
mean(attr(testResult[[1]], which="match.length", exact=TRUE)) 
min(attr(testResult[[1]], which="match.length", exact=TRUE)) 


# how would number of hits depend on number of GC regions and gaps between
gaps <- seq(0,30, by=5)
n.hits <- seq(10,200, by=5)

count.matrix <- matrix(0, nrow=length(gaps), ncol=length(n.hits), dimnames=list(gaps=as.character(gaps), n.hits=paste0("X",n.hits)))
med.matrix <- matrix(0, nrow=length(gaps), ncol=length(n.hits), dimnames=list(gaps=as.character(gaps), n.hits=paste0("X",n.hits)))
mean.matrix <- matrix(0, nrow=length(gaps), ncol=length(n.hits), dimnames=list(gaps=as.character(gaps), n.hits=paste0("X",n.hits))) 
max.matrix <- matrix(0, nrow=length(gaps), ncol=length(n.hits), dimnames=list(gaps=as.character(gaps), n.hits=paste0("X",n.hits)))
min.matrix <- matrix(0, nrow=length(gaps), ncol=length(n.hits), dimnames=list(gaps=as.character(gaps), n.hits=paste0("X",n.hits)))

for(i in 1:length(gaps))  {
  thisGap <- gaps[i]
  print(i)
  
  for(j in 1:length(n.hits))  {
	cat(j)
      cat(" . ")
    thisHits <- n.hits[j]
    pattern <- paste0("(CG.{0,", thisGap, "}){", thisHits, "}CG")
    grepResult <- gregexpr(pattern, thisChrom)
    count.matrix[i,j] <- length(grepResult[[1]])
	med.matrix[i,j] <- median(attr(grepResult[[1]], which="match.length", exact=TRUE))
	mean.matrix[i,j] <- mean(attr(grepResult[[1]], which="match.length", exact=TRUE))
	min.matrix[i,j] <- min(attr(grepResult[[1]], which="match.length", exact=TRUE))
	max.matrix[i,j] <- max(attr(grepResult[[1]], which="match.length", exact=TRUE))


  }
}

write.table(count.matrix, file="data/metrics_pc/CpG_island_size_counts.tab", sep="\t", quote=F, col.names=NA)
write.table(med.matrix, file="data/metrics_pc/CpG_island_size_median.tab", sep="\t", quote=F, col.names=NA)
write.table(mean.matrix, file="data/metrics_pc/CpG_island_size_mean.tab", sep="\t", quote=F, col.names=NA)
write.table(max.matrix, file="data/metrics_pc/CpG_island_size_max.tab", sep="\t", quote=F, col.names=NA)
write.table(min.matrix, file="data/metrics_pc/CpG_island_size_min.tab", sep="\t", quote=F, col.names=NA)

# calculate the maximum and minimum possible match lengths.
maxPoss.matrix <-  matrix(0, nrow=length(gaps), ncol=length(n.hits), dimnames=list(gaps=as.character(gaps), n.hits=paste0("X",n.hits)))
minPoss.matrix <-  matrix(0, nrow=length(gaps), ncol=length(n.hits), dimnames=list(gaps=as.character(gaps), n.hits=paste0("X",n.hits)))
for(i in 1:length(gaps))  {
  thisGap <- gaps[i]
  print(i)
  
  for(j in 1:length(n.hits))  {
	thisHits <- n.hits[j]
	minPoss.matrix[i,j] <- (thisHits  * 2)  + 2
	maxPoss.matrix[i,j] <- (thisHits  * (2+thisGap)) + 2

	
  }
}
write.table(maxPoss.matrix, file="data/metrics_pc/CpG_island_size_maxPoss.tab", sep="\t", quote=F, col.names=NA)
write.table(minPoss.matrix, file="data/metrics_pc/CpG_island_size_minPoss.tab", sep="\t", quote=F, col.names=NA)


#count.matrix <- read.delim("data/metrics_pc/CpG_island_size_counts.tab")

plot(count.matrix[,1]~ gaps, type="l", ylab="count", xlab="gap length")
for(i in 2:ncol(count.matrix)) {
  lines(count.matrix[,i] ~ gaps)
}
for(i in 3:8)  {
  text(75, y=count.matrix["75", paste0("X",i)], labels = paste0("x", i), pos = 2)
}
text(75, y=count.matrix["75", "X4"], labels = "x4", pos = 2)
text(75, y=count.matrix["75", "X5"], labels = "x5", pos = 2)



# how do the median values compare to the maxPoss values - most are <0.6
 med.matrix / maxPoss.matrix
plot(med.matrix[,1] / maxPoss.matrix[,1]~ gaps, type="l", ylab="count", xlab="gap length", ylim=c(0,1))
for(i in 2:ncol(count.matrix)) {
  lines(med.matrix[,i] / maxPoss.matrix[,i] ~ gaps)
}


# plot max hit lenght as ratio of max possible 
png(file="data/metrics_pc/CpG_island_size_maxOverMaxPoss.png")
plot(max.matrix[,1] / maxPoss.matrix[,1]~ gaps, type="l", ylab="count", xlab="gap length", ylim=c(0,1), main="Max / MaxPoss length")
for(i in 2:ncol(count.matrix)) {
  lines(max.matrix[,i] / maxPoss.matrix[,i] ~ gaps)
}
for(i in c(40,50,60,70,80,90,100))  {
  text(15, y=max.matrix["15", paste0("X",i)] / maxPoss.matrix["15", paste0("X",i)] , labels = paste0("x", i), pos =4, cex=.7)
}
dev.off()




# Gardiner-Garden and XX 1987 wrote a definition for CpG islands  O/E > 0.6 and GC% > 50  (over 100bp)
# Later expanded and used by others at 200bp or 500bp
# Seem's fairly arbitrary but may be useful  to compare with other methods.
# Their formula for O/E  was   N.CpG  x   L   /  (N.C  x N.G)    

#So, for a 200bp stretch with 51 G and 51 C (GC% > 50%), 
 N.CpG <- (.6 /200)  * (51 * 51)    # only 7.803

# with GC% at 60%, 
 N.CpG <- (.6 /200)  * (60 * 60)  # 10.8

countPattern("CG", testString)

letterFrequency(DNAString(testString), letters = c("C","G"))

nchar(testString)

letterFrequencyInSlidingView(DNAString(testString), letters = c("C","G"), view.width = 200)



cpgRatio <- function(x, verbose=FALSE)  {
  # calculate the O/E CpG ratio for a give sequence
  # would be good to be able to apply this in a sliding window.
  x <- DNAString(x)
  L <- length(x)
  N.CpG <- countPattern("CG", x)
  lf.CG <- letterFrequency(x, letters = c("C","G"))
  N.C <- lf.CG[1]
  N.G <- lf.CG[2]
  GCpc <- sum(lf.CG/L) * 100
  if(verbose)  {
    cat(paste("Length:", L, "\n", "N.CpG", N.CpG, "\n", "N.C", N.C, "N.G", N.G, "\n", "GC%", format(GCpc, digits=3), "%", "\n"))
  }
  return(as.numeric( (N.CpG  *   L)   /  (N.C  * N.G)))
  
}

cpgRatio(testString, verbose=T)





results3 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland)  # from mus.x-inactivationModel.R
# how many of these have CpG ratio above .6 or GC% above 50%?

# many are not 200bp, but they are delimited by CpG so should be good
seqs <- getSeq( genome, GRanges("chrX", results3))     # super cool that this works!
cpgRatioVec <- numeric()
cgPcVec <- numeric()
for(i in 1:length(results3)) {
  cpgRatioVec[i] <- cpgRatio(seqs[[i]])
  cgPcVec  [i] <- sum(letterFrequency(seqs[[i]], letters = c("C","G")))/length(seqs[[i]])
  
}

table(cpgRatioVec > 0.6, cgPcVec > 0.5)  # bf.CpGisland is VERY good at finding CpG islands.

# the only problem may be the length

hist(width(results3))
hist(cpgRatioVec, breaks=100)
hist(cgPcVec, breaks=100)  # odd spike at ~75%   are these repeats?
