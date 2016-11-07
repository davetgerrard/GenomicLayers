#mus.plotSimonXinactivation.R
setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Mmusculus.UCSC.mm9) 
genome <- BSgenome.Mmusculus.UCSC.mm9 
# digression from mus.x-inactivationModel.R  - compare with real mouse data

OUTPUTDIR <- "results/mus_X_inactivation/"
GENOME <- 'mm9'


# read in bedGraph Data
fileName <- "C:/Users/Dave/data/Simon_etal_2013_MusXist/bedGraphs/GSM1182890_d3.xist.mus.bedGraph.gz"

bg <- read.table(fileName, header=F, skip=1, sep= " ")

nrow(bg)
names(bg) <- c("chrom", "start", "end", "value")
# only want chrX values.
bg.chrom <- subset(bg , chrom=="chrX")
nrow(bg.chrom)


  
  
# plot it !
library(Sushi)

chrom   <-"chrX"
chromstart <- 1
chromend <-  length(genome[["chrX"]] )
#plotBedgraph(bg.chrom, chrom="chrX", chromstart=1, chromend=170000000)  # works

png(paste0(OUTPUTDIR, "x.inactivation.simonData.png"), width=2400, height=700)
plotBedgraph(bg.chrom, chrom="chrX", chromstart=chromstart, chromend=chromend)  # works
#plotBedgraph(bg.chrom,chrom,chromstart,chromend,colorbycol= SushiColors(5))  #does not work
labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
dev.off()


# measure correlation with simulation in different waves. 


#load the simulation results
n.waves <- 1000
n.waves <- 10
n.iters <- 1000
load(paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves,  ".", n.iters,".waveList.Rdata"))  # may not work for some earlier runs

# simon et al., data are in 500bp bins. 
# Could get measure of how much of each bin is marked in simulation waves, then do correlation

# simplest would be whether regions are overlapped or not. 0/1


# need to take a little care of the chromosome start and end, which appear non-binned in the Simon et al. data

nrow(bg.chrom)   # many rows.
bg.chrom.GR <- GRanges(bg.chrom$chrom, IRanges(start=bg.chrom$start, end=bg.chrom$end), value=bg.chrom$value)
bg.chrom.IR <- IRanges(start=bg.chrom$start, end=bg.chrom$end)

thisFactor <- "H3K27me3"

ov <- overlapsAny(bg.chrom.IR, waveList[[10]][['layerSet']][[thisFactor]])
sum(ov)
length(ov)

boxplot(as.numeric(bg.chrom$value) ~ ov)
by(as.numeric(bg.chrom$value), ov, median)

cor(bg.chrom$value, as.integer(ov))

cor.wl <- numeric()
cor.wl[1] <- 0
for(i in 2:n.waves)  {
  ov <- overlapsAny(bg.chrom.IR, waveList[[i]][['layerSet']][[thisFactor]])
  cor.wl[i] <- cor(bg.chrom$value, as.integer(ov))
}
png(paste0(OUTPUTDIR, "x.inactivation.simonData.corByWave.", n.waves, ".", n.iters,".png"), width=2400, height=700, res=150)
plot(1:length(cor.wl),cor.wl  , xlab="wave", ylab="correlation", main=thisFactor)
dev.off()
# not great correlations, but does seem to changes systmatically.


# count how many blocks above a threshold are hit in each wave.

threshold <- .5

bg.chrom.subset <- subset(bg.chrom, value >=threshold)
bg.chrom.subset.IR <- IRanges(start=bg.chrom.subset$start, end=bg.chrom.subset$end)
n.subset <- nrow(bg.chrom.subset)
bg.chrom.inverse <- subset(bg.chrom, value <threshold)
bg.chrom.inverset.IR <- IRanges(start=bg.chrom.inverse$start, end=bg.chrom.inverse$end)

count.wl <- integer()
count.inverse <- integer()   # also count how many below threshold are hit.
count.wl[1] <- count.inverse [1] <- 0
for(i in 2:n.waves)  {
  ov <- overlapsAny(bg.chrom.subset.IR, waveList[[i]][['layerSet']][[thisFactor]])
  count.wl[i] <- sum(ov)
  ov <- overlapsAny(bg.chrom.inverset.IR, waveList[[i]][['layerSet']][[thisFactor]])
  count.inverse[i] <- sum(ov)
}

png(paste0(OUTPUTDIR, "x.inactivation.simonData.hitCounts.", n.waves, ".", n.iters,".png"), width=1200, height=1200, res=150)
par(mfrow=c(3,1))
plot(1:length(count.wl),count.wl  , xlab="wave", ylab="hit count", main=thisFactor)
plot(1:length(count.inverse),count.inverse  , xlab="wave", ylab="miss count", main=thisFactor)
hit.miss <- count.wl - count.inverse
plot(1:length(hit.miss),hit.miss  , xlab="wave", ylab="hit - miss", main=thisFactor)
dev.off()

# nope, that first run really sucks!  - only the very early hits (e.g. CpG_islands come close to the real data)




# try using different measure, as per TSS optimisation methods
test_function <- function(layerList, targetLayer=target.layer, target.vec)  {
  inter.size <- sum(width(intersect(layerList$layerSet[[targetLayer]], target.vec)))
  union.size <- sum(width(union(layerList$layerSet[[targetLayer]], target.vec)))
  #layer.vec <- as.numeric(strsplit(as.character(layerList$layerSet[[targetLayer]]),"")[[1]])
  return(inter.size/ union.size)
}


threshold <-.5

bg.chrom.subset <- subset(bg.chrom, value >=threshold)
bg.chrom.subset.IR <- IRanges(start=bg.chrom.subset$start, end=bg.chrom.subset$end)
n.subset <- nrow(bg.chrom.subset)
bg.chrom.inverse <- subset(bg.chrom, value <threshold)
bg.chrom.inverset.IR <- IRanges(start=bg.chrom.inverse$start, end=bg.chrom.inverse$end)



test_function(waveList[[10]], targetLayer=thisFactor, target.vec=bg.chrom.subset.IR)


score.wl <- integer()
score.inverse <- integer()   # also count how many below threshold are hit.
score.wl[1] <- score.inverse [1] <- 0
for(i in 2:n.waves)  {

  score.wl[i] <- test_function(waveList[[i]], targetLayer=thisFactor, target.vec=bg.chrom.subset.IR)

  score.inverse[i] <- test_function(waveList[[i]], targetLayer=thisFactor, target.vec=bg.chrom.inverset.IR)
}


png(paste0(OUTPUTDIR, "x.inactivation.simonData.testScores.", n.waves, ".", n.iters,".png"), width=1200, height=1200, res=150)
par(mfrow=c(3,1))
plot(1:length(score.wl),score.wl  , xlab="wave", ylab="score", main=thisFactor)
plot(1:length(score.inverse),score.inverse  , xlab="wave", ylab="miss score", main=thisFactor)
hit.miss <- score.wl - score.inverse
plot(1:length(hit.miss),hit.miss  , xlab="wave", ylab="hit - miss", main=thisFactor)
dev.off()



test_function(waveList[[i]], targetLayer=thisFactor, target.vec=waveList[[i]][['layerSet']][['CpG_island']])
#only 5%

test_function(waveList[[i]], targetLayer="CpG_island", target.vec=bg.chrom.subset.IR)

test_function(waveList[[5]], targetLayer="CpG_island", target.vec=bg.chrom.subset.IR)

test_function(waveList[[i]], targetLayer="PRC", target.vec=bg.chrom.subset.IR)

test_function(waveList[[5]], targetLayer="PRC", target.vec=bg.chrom.subset.IR)


#  Poor agreement between the simulation and the Simon et al data. 
# I wonder if the simulated CpG_islands may be a better fit OR whether they match the simulated H3K27me3 (wondering whether the offset parameter is too generours).

window.start <- 60000000
window.end <- 140000000

window.start <- 1
window.end <- length(thisChrom)


layerSubset <- restrict(waveList[[10]][['layerSet']][[thisFactor]], start = window.start, end=window.end)
bgSubset <- restrict(bg.chrom.GR, start = window.start, end=window.end)
cpgHitsSubset <- restrict(results3, start = window.start, end=window.end)

png(filename=paste0(OUTPUTDIR, "x.inactivation.simonData.zoomWindow.", n.waves, ".", n.iters,".", ceiling(window.start/1000000), ".", ceiling(window.end/1000000),".png"), width=2400, height=800)
#par(mfrow=c(2,1))
#plot(start(restrict(waveList[[10]][['layerSet']][[thisFactor]], start = window.start, end=window.end)), start(restrict(waveList[[10]][['layerSet']][[thisFactor]], start = window.start, end=window.end)) )
plot(start(bgSubset), bgSubset$value)
grid(nx=40, ny=0, col="grey50")
#points(start(layerSubset), rep(10, length(layerSubset)))
rug(start(layerSubset), side=3, line=2)
mtext("final layer", side=3, line=1, col="black", adj = 0)
rug(start(cpgHitsSubset), side=3, line=4, col="blue")
mtext("CpG islands", side=3, line=3, col="blue", adj=0)
# this Xist regions is highly punctate.    using 
#window.start <- 100000000
#window.end <- 105000000
dev.off()
# using 80Mb to 120Mb,  there is a high density of CpG hits (
length(results3)  # 2148
length(cpgHitsSubset) # 489


# same plot showing more waves in wider upper margin

wavesToShow <- c(2,4,6,8,10)
window.start <- 1
window.end <- length(thisChrom)

png(filename=paste0(OUTPUTDIR, "x.inactivation.simonData.zoomWindow.", n.waves, ".", n.iters,".", ceiling(window.start/1000000), ".", ceiling(window.end/1000000),".wavesRug.png"), width=2400, height=800)
par(mar=c(5,5,(length(wavesToShow) + 2)*2,4))
plot(start(bgSubset), bgSubset$value, ylab="Xist signal", xlab="mm9 chrX")
grid(nx=40, ny=0, col="grey50")
for(i in 1:length(wavesToShow)) {
  layerSubset <- restrict(waveList[[wavesToShow[i]]][['layerSet']][[thisFactor]], start = window.start, end=window.end)
  plotLine <- (i + 1)*2
  rug(start(layerSubset), side=3, line=plotLine)
  mtext(paste("Wave", wavesToShow[i]), side=3, line=plotLine, col="black", adj = 0)
}
rug(start(cpgHitsSubset), side=3, line=2, col="blue")
mtext("CpG islands", side=3, line=2, col="blue", adj=0)
dev.off()






# TODO - try concordance as a measure?
