#mus.plotSimonXinactivation.R


# digression from mus.x-inactivationModel.R  - compare with real mouse data

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


# simon et al., data are in 500bp bins. 
# Could get measure of how much of each bin is marked in simulation waves, then do correlation

# simplest would be whether regions are overlapped or not. 0/1


# need to take a little care of the chromosome start and end, which appear non-binned in the Simon et al. data

nrow(bg.chrom)   # many rows.
bg.chrom.GR <- GRanges(bg.chrom$chrom, IRanges(start=bg.chrom$start, end=bg.chrom$end))
bg.chrom.IR <- IRanges(start=bg.chrom$start, end=bg.chrom$end)

ov <- overlapsAny(bg.chrom.IR, waveList[[100]][['layerSet']][[thisFactor]])
sum(ov)
length(ov)

boxplot(as.numeric(bg.chrom$value) ~ ov)

cor(bg.chrom$value, as.integer(ov))

cor.wl <- numeric()
cor.wl[1] <- 0
for(i in 2:200)  {
  ov <- overlapsAny(bg.chrom.IR, waveList[[i]][['layerSet']][[thisFactor]])
  cor.wl[i] <- cor(bg.chrom$value, as.integer(ov))
}
png(paste0(OUTPUTDIR, "x.inactivation.simonData.corByWave.png"), width=2400, height=700, res=150)
plot(1:length(cor.wl),cor.wl  , xlab="wave", ylab="correlation", main=thisFactor)
dev.off()
# not great correlations, but does seem to changes systmatically.


# count how many blocks above a threshold are hit in each wave.

threshold <- .25

bg.chrom.subset <- subset(bg.chrom, value >=threshold)
bg.chrom.subset.IR <- IRanges(start=bg.chrom.subset$start, end=bg.chrom.subset$end)
n.subset <- nrow(bg.chrom.subset)
bg.chrom.inverse <- subset(bg.chrom, value <threshold)
bg.chrom.inverset.IR <- IRanges(start=bg.chrom.inverse$start, end=bg.chrom.inverse$end)

count.wl <- integer()
count.inverse <- integer()   # also count how many below threshold are hit.
count.wl[1] <- count.inverse [1] <- 0
for(i in 2:200)  {
  ov <- overlapsAny(bg.chrom.subset.IR, waveList[[i]][['layerSet']][[thisFactor]])
  count.wl[i] <- sum(ov)
  ov <- overlapsAny(bg.chrom.inverset.IR, waveList[[i]][['layerSet']][[thisFactor]])
  count.inverse[i] <- sum(ov)
}

png(paste0(OUTPUTDIR, "x.inactivation.simonData.hitCounts.png"), width=1200, height=1200, res=150)
par(mfrow=c(3,1))
plot(1:length(count.wl),count.wl  , xlab="wave", ylab="hit count", main=thisFactor)
plot(1:length(count.inverse),count.inverse  , xlab="wave", ylab="miss count", main=thisFactor)
hit.miss <- count.wl - count.inverse
plot(1:length(hit.miss),hit.miss  , xlab="wave", ylab="hit - miss", main=thisFactor)
dev.off()

# nope, that first run really sucks!  - only the very early hits (e.g. CpG_islands come close to the real data)














