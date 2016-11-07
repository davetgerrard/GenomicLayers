
# Simulate 2-stage X-inactivation in mouse
#environment, packages and functions

# Rscript.exe mus.x-inactivationModel.fromConfig.R --config CONFIG-FILE.R
.libPaths(c(.libPaths(), "C:/Users/Dave/Documents/R/win-library/3.3"))  # required to run from command line and pick up local package library
# Load CONFIG --------------------------------
library(getopt)
### this file has been setup to use a config script specifying the input and output files. 

spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'config'  , 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec);
if ( is.null(opt$config ) ) { opt$config <- 'C:/Users/Dave/HalfStarted/predictFromSequence/scripts/CONFIG.mus.x-inactivationModel.2016-10-27.R'}


source(opt$config)

if(!file.exists(OUTPUTDIR))  dir.create(OUTPUTDIR)



#mus.plotSimonXinactivation.R
#setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#library(BSgenome.Mmusculus.UCSC.mm10)
#library(BSgenome.Mmusculus.UCSC.mm9) 
#genome <- BSgenome.Mmusculus.UCSC.mm9 
# digression from mus.x-inactivationModel.R  - compare with real mouse data

#OUTPUTDIR <- "results/mus_X_inactivation/"
#GENOME <- 'mm9'


# read in bedGraph Data
fileName <- "C:/Users/Dave/data/Simon_etal_2013_MusXist/bedGraphs/GSM1182890_d3.xist.mus.bedGraph.gz"

bg <- read.table(fileName, header=F, skip=1, sep= " ")

nrow(bg)
names(bg) <- c("chrom", "start", "end", "value")
# only want chrX values.
bg.chrom <- subset(bg , chrom=="chrX")
nrow(bg.chrom)




#load the simulation results

#load(paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves,  ".", n.iters,".waveList.Rdata"))  # may not work for some earlier runs

load( file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", further.waves, ".", n.iters,".waveList.Rdata"))


# want to plot the CpG island factor binding sites
results3 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland)



nrow(bg.chrom)   # many rows.
bg.chrom.GR <- GRanges(bg.chrom$chrom, IRanges(start=bg.chrom$start, end=bg.chrom$end), value=bg.chrom$value)
bg.chrom.IR <- IRanges(start=bg.chrom$start, end=bg.chrom$end)

#thisFactor <- "H3K27me3"

#window.start <- 1
#window.end <- length(thisChrom)


#layerSubset <- restrict(waveList[[10]][['layerSet']][[thisFactor]], start = window.start, end=window.end)
bgSubset <- restrict(bg.chrom.GR, start = window.start, end=window.end)
cpgHitsSubset <- restrict(results3, start = window.start, end=window.end)






# same plot showing more waves in wider upper margin

#wavesToShow <- c(2,4,6,8,10)


png(filename=paste0(OUTPUTDIR, "x.inactivation.simonData.zoomWindow.", further.waves, ".", n.iters,".", ceiling(window.start/1000000), ".", ceiling(window.end/1000000),".wavesRug.png"), width=2400, height=800)
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


windowWidth <- function(x, starts=NULL, window.size=NULL,  n.windows=length(starts), limit=max(starts)+1) {
  require(IRanges)
  # return the width of all features in in consecutive windows.
  # x should be an IRanges object
  
  
  
  width.vec <- integer()
  
  for(i in 1:n.windows)  {
    thisWindow.IR <- restrict(IRanges(start=starts[i], width=window.size), start = 1, end=limit)
    width.vec[i] <- sum(width(intersect(x, thisWindow.IR)))
    
  }
  return(width.vec)
  
}

window.length <- 1*1000*1000

bin.starts <- seq(1, window.end, by=0.1*1000*1000)

#widthVec3 <- windowWidth(results3, starts=bin.starts, window.size=window.length, limit=window.end)
#widthVecWave <- windowWidth(waveList[[200]][['layerSet']][[thisFactor]], starts=bin.starts, window.size=window.length, limit=window.end)
#plot(1:length(widthVec3),widthVec3, type="l")

png(filename=paste0(OUTPUTDIR, "x.inactivation.simonData.zoomWindow.", further.waves, ".", n.iters,".", ceiling(window.start/1000000), ".", ceiling(window.end/1000000),".sliding.png"), width=2400, height=800)
par(mfrow=c(2,1))
layerSubset <- restrict(waveList[[max(wavesToShow)]][['layerSet']][[thisFactor]], start = window.start, end=window.end)
widthVecWave <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
plot(1:length(widthVecWave),widthVecWave, type="l")
for(i in 1:length(wavesToShow)) {
  layerSubset <- restrict(waveList[[wavesToShow[i]]][['layerSet']][[thisFactor]], start = window.start, end=window.end)
  widthVecWave <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
  lines(1:length(widthVecWave),widthVecWave)
  #plotLine <- (i + 1)*2
  #rug(start(layerSubset), side=3, line=plotLine)
  #mtext(paste("Wave", wavesToShow[i]), side=3, line=plotLine, col="black", adj = 0)
}
plot(start(bgSubset), bgSubset$value, ylab="Xist signal", xlab="mm9 chrX")
dev.off()

# TODO - fix the above so that they are properly aligned!  or separate into two.

#  TODO - calculate mean bedgraph score for the mus data using modified version of windowWidth
#           then can properly correlate scores for each wave. 


windowMean <- function(x, y,  starts=NULL, window.size=NULL,  n.windows=length(starts), 
                         window.func = mean, limit=max(starts)+1) {
  require(IRanges)
  # return the value of function on a certain score column in in consecutive windows.
  # x should be an IRanges object
  # y should be a vector of values of same lenght as x
  #     for several reasons, cannot yet use mcols() on x
  
  # there should be no self overlaps (otherwiise cannot make fair averages). 
  
  stopifnot( length(findOverlaps(x, drop.self=TRUE) ) ==0)   # should be no non-self overlaps.
  stopifnot(length(x) == length(y))
  
  mean.vec <- numeric()
  
  for(i in 1:n.windows)  {
    thisWindow.IR <- IRanges(start=starts[i], width=window.size)
    window.index <- overlapsAny(x, thisWindow.IR)
    
    # these two not currently used, was thinking to weight by number of regions or amount of space taken by regions within a window.
    # Will matter for uneven sized regions.
    window.coverage <- sum(width(x[window.index]))
    window.count <- sum(window.index)
    
    mean.vec[i] <- mean(y[window.index])
    
  }
  return(mean.vec)
  
}


meanScores <- windowMean(bg.chrom.IR, y=mcols(bg.chrom.GR)$value, starts=bin.starts, window.size=window.length, limit=window.end)
hist(meanScores)

length(meanScores)
length(widthVecWave)

cor(meanScores, widthVecWave)


plot(1:length(widthVecWave),widthVecWave, type="l")
plot(1:length(meanScores),meanScores, type="l")


# could now plot these on same axis - but need to scale all wave values by the same max.
# so,  calculate all waveVecs first, then go through plotting second.

waveVecList <- list()
cor.vec <- numeric()
max.value <- NULL
no.index <- seq(30, length(meanScores), by=10)  # evenly spaced, non-overlapping windows  
for(i in 1:length(wavesToShow)) {
  layerSubset <- restrict(waveList[[wavesToShow[i]]][['layerSet']][[thisFactor]], start = window.start, end=window.end)
  waveVecList[[i]] <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
  max.value <- max(max.value, max(waveVecList[[i]]))
  cor.vec[i] <- cor(meanScores[no.index], waveVecList[[i]][no.index] , use="complete.obs", method="spearman")

}



png(filename=paste0(OUTPUTDIR, "x.inactivation.simonData.zoomWindow.", further.waves, ".", n.iters,".", ceiling(window.start/1000000), ".", ceiling(window.end/1000000),".sliding.MeanBins.png"), width=2400, height=1000, res=150)
#par(mfrow=c(2,1))
nonNaVec <- which(!is.na(meanScores))
plot(nonNaVec,meanScores[nonNaVec]/max(meanScores[nonNaVec]), type="l", col="blue",lwd=2, xlim=c(0, length(meanScores)), xlab="100kb bin", ylab="signal")  # max scaled
#layerSubset <- restrict(waveList[[max(wavesToShow)]][['layerSet']][[thisFactor]], start = window.start, end=window.end)
#widthVecWave <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
#plot(1:length(widthVecWave),widthVecWave, type="l")
for(i in 1:length(wavesToShow)) {
  #layerSubset <- restrict(waveList[[wavesToShow[i]]][['layerSet']][[thisFactor]], start = window.start, end=window.end)
  widthVecWave <- waveVecList[[i]]
  
  greyColour <- paste0("gray",round(60 - ((wavesToShow[i]/max(wavesToShow))*55)))
  lines(1:length(widthVecWave),widthVecWave/max.value, col=greyColour)
  #plotLine <- (i + 1)*2
  #rug(start(layerSubset), side=3, line=plotLine)
  #mtext(paste("Wave", wavesToShow[i]), side=3, line=plotLine, col="black", adj = 0)
}
# replot the meanScores to  make it more clear
#lines(nonNaVec,meanScores[nonNaVec]/max(meanScores[nonNaVec]), col="blue",lwd=2)  # max scaled
#plot(start(bgSubset), bgSubset$value, ylab="Xist signal", xlab="mm9 chrX")
dev.off()



# DONE - is correlation correct when using overlapping sliding windows?
# perhaps need to measure using non-overlapping subset.
# Currently using 1MBp windows every 100kb.

cor(meanScores, waveVecList[[i]] , use="complete.obs", method="spearman")
no.index <- seq(30, length(meanScores), by=10)  # non-overlapping windows
cor(meanScores[no.index], waveVecList[[i]][no.index] , use="complete.obs", method="spearman")

plot(wavesToShow, cor.vec, ylim=c(-1,1))
plot(wavesToShow, cor.vec)




