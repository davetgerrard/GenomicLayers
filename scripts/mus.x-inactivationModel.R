

# Simulate 2-stage X-inactivation in mouse
#environment, packages and functions


# TODO - CpG BF finding CpGs, but results not very life-like. Check paper, is it GC% instead?
# TODO - better/best method to assess results?

setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Mmusculus.UCSC.mm9)  # to compare with data from Simon et al.
#require(Biostrings)   # included in mm10 package dependencies
source('scripts/pfs.functions.R')

OUTPUTDIR <- "results/mus_X_inactivation/"
GENOME <- 'mm9'

# load X chromosome sequence 
genome <- BSgenome.Mmusculus.UCSC.mm9 
thisChrom <- genome[["chrX"]] 

# set up a layerSet on chrom X.

layerSet.X <- list(LAYER.0 = thisChrom)
layerSet.X[['CpG_island']] <- IRanges()
layerSet.X[['PRC']] <- IRanges()
layerSet.X[['H3K27me3']] <- IRanges()

layerList.X <- list(layerSet=layerSet.X, history=NULL)  # add some metadata



# specify binding factors in model (and proportions?)

#   Stage 1.  - Gene rich regions  
#               CpG islands    [Pinter et al., 2012](http://www.citeulike.org/user/daveGerrard/article/11359142)
#   Stage 2. - Spread out from these regions in a two stage feedback (could start with one factor?).
#           - Meant to represent recruitment of PRC and then marking of H3K27me3.


bf.CpG <- createBindingFactor.DNA_motif("CpG", patternString="CG", 
                                        profile.layers=NULL,profile.marks=NULL, 
                                        mod.layers = "CpG_island", mod.marks=1)

bf.PRC <- createBindingFactor.layer_region("PRC",  patternLength=200, mismatch.rate=0, 
                                           profile.layers = "CpG_island", profile.marks = 1,  
                                           mod.layers = "PRC", mod.marks=1, stateWidth=500)

#bf.CpGisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){50}CG", patternLength=300, profile.layers="H3K27me3",profile.marks=0, mod.layers = "CpG_island", mod.marks=1)
# needed more seed areas, so took shorter CpG motifs.
# bf.CpGisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){9}CG", patternLength=20,
#                                                profile.layers="H3K27me3",profile.marks=0, 
#                                                mod.layers = "CpG_island", mod.marks=1, stateWidth=300)
bf.CpGisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){9}CG", patternLength=20,
                                               mod.layers = "CpG_island", mod.marks=1, stateWidth=200)

#N.B. temporarily setting a fixed patternLength, because don't know how to implement variable patternLength yet.
# current effect is that hits < patternLength are discarded

# want a binding factor that creates mods either up or downstream of binding site. New parameter 'offset.method' to accept function as argument
upDownFunc <- function(x)  {
  return(sample(c(x, -x), 1))
  
}

bf.spreadRep <- createBindingFactor.layer_region("spreadRep", patternLength=150, mismatch.rate=0, 
                                                 profile.layers = "PRC", profile.marks = 1, 
                                                 mod.layers = "H3K27me3", mod.marks=1, 
                                                 stateWidth=200,offset=350, offset.method=upDownFunc)
# 2016-10-19 reduced the statewidth and offset because I thought it was growing too fast relative to CpG islands. 


results1 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpG)    # should be many hits, each of 2bp.  
results2 <- matchBindingFactor(layerSet.X, bindingFactor = bf.PRC)    # should be no hits.
results3 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland)
results4 <- matchBindingFactor(layerSet.X, bindingFactor = bf.spreadRep)    # should be no hits.

results4b <- matchBindingFactor(layerList.X$layerSet, bindingFactor = bf.spreadRep)  



##!!!!
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
#XFS <- list(bf.CpGisland =bf.CpGisland , bf.PRC.1=bf.PRC, bf.spreadRep.1=bf.spreadRep, bf.PRC.2=bf.PRC, bf.spreadRep.2=bf.spreadRep)
XFS <- list(bf.CpGisland =bf.CpGisland , bf.PRC.1=bf.PRC, bf.spreadRep.1=bf.spreadRep, bf.PRC.2=bf.PRC, bf.spreadRep.2=bf.spreadRep)

mod.X <- runLayerBinding(layerList.X, factorSet = XFS, iterations = 10000, verbose=T)  # will saturate CpG islands in first round

# these were some tests I was doing when binding wasn't working properly. Changed matchBindingFactor to use intersect()
#mod.CpGIsland <- runLayerBinding(layerList.X, factorSet = list(bf.CpGisland =bf.CpGisland), iterations = 10000, verbose=T)
#mod.PRC <- runLayerBinding(mod.CpGIsland, factorSet = list(bf.PRC=bf.PRC), iterations = 10000, verbose=T)   # SHOULD be only a small number in relation to CpG islands.
#resultsPRC <- matchBindingFactor(mod.CpGIsland$layerSet, bindingFactor = bf.PRC) 


# TODO would be good to add turning off of CpG islands as repression spreads. Representing methylation. Could be modelled as separate layer or 




# run several waves of XFS to show spreading.  
# need to collect data each time to plot. it. 

n.waves <- 10
#n.waves <- 100
n.iters <- 2000
waveList <- list()
waveList[[1]] <- current.X <- layerList.X
for(i in 2:n.waves) {
  print(paste("Running wave", i))
  current.X <- runLayerBinding(current.X, factorSet = XFS, iterations = n.iters, verbose=T)
  waveList[[i]] <- current.X
}


save(XFS, waveList, file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves, ", ", n.iters,".waveList.Rdata"))

# 
## run further waves - uncomment to run, takes a fair while. Saves every 100 waves.
# n.waves <- 100
# load(paste0(OUTPUTDIR, "x.inactivation.", n.waves, ".waveList.Rdata"))
current.X <- waveList[[n.waves]]
next.wave <- n.waves + 1
further.waves <- 1000
for(i in next.wave:further.waves) {
  print(paste("Running wave", i))
  current.X <- runLayerBinding(current.X, factorSet = XFS, iterations = 100, verbose=T)
  waveList[[i]] <- current.X
  if(i %% 100 == 0)  {
    save(XFS, waveList, file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", i, ".waveList.Rdata"))
  }
}

save(XFS, waveList, file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", further.waves, ".waveList.Rdata"))



width.waves <- integer()
for(i in 1:n.waves) {
  width.waves[i] <- sum(width(waveList[[i]]$layerSet$H3K27me3))
}

png(paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves, ".H3K27me3.png"))
plot(1:n.waves, width.waves)  # why is this exact?   (because the number of domains is limited to same number ecah time?
dev.off()

#plot(start(waveList[[i]]$layerSet$CpG_island),start(waveList[[i]]$layerSet$H3K27me3))

i <- 2
i <- n.waves
par(mfrow=c(2,1))
hist(start(waveList[[i]]$layerSet$CpG_island), breaks=500)
hist(start(waveList[[i]]$layerSet$H3K27me3), breaks=500)
# compare with x-inactivation data.


thisFactor <- "CpG_island"
for(thisFactor in c( "CpG_island" ,"PRC"   ,"H3K27me3")) {
  pdf(paste0(OUTPUTDIR, "x.inactivation.", GENOME,".",n.waves,".", thisFactor,".each10.pdf"), height=20)
  par(mfrow=c(11,1))
  for(k in seq(10, 100, by=10)) {
    hist(start(waveList[[k]][['layerSet']][[thisFactor]]), breaks=500, main=paste(thisFactor, "wave", k), xlab="X chromosome")
  }
  dev.off()
}
#optimise binding factors to better match sequential inactivation

# output further waves if available  - very little difference between runs 100 and 1000 (of first attempt)
for(thisFactor in c( "CpG_island" ,"PRC"   ,"H3K27me3")) {
  pdf(paste0(OUTPUTDIR, "x.inactivation.", GENOME,".",further.waves,".", thisFactor,".each10.pdf"), height=20)
  par(mfrow=c(10,1))
  for(k in seq(100, 1000, by=100)) {
    hist(start(waveList[[k]][['layerSet']][[thisFactor]]), breaks=500, main=paste(thisFactor, "wave", k), xlab="X chromosome")
  }
  dev.off()
}


# for a quick binning method with ggplot2 https://www.biostars.org/p/69748/
# How about a chromosome wide heatmap plus genes?




library(IdeoViz)
require(RColorBrewer)

#from IdeoViz vignette

ideo_mm <- getIdeo(GENOME)
#chroms <- c("chr1","chr2","chrX")
thisChrom <- "chrX"

chrom_bins <- getBins(thisChrom, ideo_mm,stepSize=1*1000*1000)
focal.waves <- c(2,4,8,16,32,64)  # too similar for H3K27me3
focal.waves <- c(10,25,50,75,100)
focal.waves <- c(2,4,6,8,10)

for(thisFactor in c( "CpG_island" ,"PRC"   ,"H3K27me3")) {

# create a data frame on chrom_bins and get coverage counts for each bin for each wave (or selected waves).
chrom.data <- as.data.frame(chrom_bins)
for(k in focal.waves)  {
  wave.start <- GRanges(thisChrom, IRanges(start=start(waveList[[k]][['layerSet']][[thisFactor]]), end=start(waveList[[k]][['layerSet']][[thisFactor]])), mcols=data.frame(score=width(waveList[[k]][['layerSet']][[thisFactor]])))
  ov <- findOverlaps(chrom_bins, wave.start)   # use start so each can only overlap one bin
  indexCounts <- tapply(mcols(wave.start)[,1], queryHits(ov), FUN=sum)
  chrom.data[,paste0("wave.", k)] <- 0
  chrom.data[names(indexCounts),paste0("wave.", k)] <- indexCounts
}

chrom.data.GR <- as(chrom.data, "GRanges")

par(mar=c(15,15,15,15))

y.max <- max(mcols(chrom.data.GR)[,ncol(mcols(chrom.data.GR))])
pdf(paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves,".", thisFactor,".plotOnIdeo.pdf"), width=13, height=7)
par(mar=c(15,15,15,15))  #not working ?
plotOnIdeo(chrom=thisChrom, # which chrom to plot?
           ideoTable=ideo_mm, # ideogram name
           values_GR=chrom.data.GR, # data goes here
           value_cols=colnames(mcols(chrom.data.GR)), # col to plot
           col=rev(brewer.pal(n= ncol(mcols(chrom.data.GR)), 'Blues')), # colours
           val_range=c(0,y.max), # set y-axis range
           ylab="Coverage",
           plot_title=thisFactor)
dev.off()
#  display.brewer.all(n=10, exact.n=FALSE)

#can zoom in using bplim parameter.
pdf(paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves,".", thisFactor,".plotOnIdeo.zoomXist.pdf"), width=13, height=7)
plotOnIdeo(chrom=thisChrom, # which chrom to plot?
           ideoTable=ideo_mm, # ideogram name
           values_GR=chrom.data.GR, # data goes here
           value_cols=colnames(mcols(chrom.data.GR)), # col to plot
           col=rev(brewer.pal(n= ncol(mcols(chrom.data.GR)), 'Blues')), # colours
           val_range=c(0,y.max), # set y-axis range
           ylab="Coverage",
           bpLim = c(80000000, 130000000),
           plot_title=thisFactor)
dev.off()

}








# DEVELOPMENT ---------------------- 



binstats <- avgByBin(featureData= as.data.frame( GRanges(thisChrom, waveList[[k]][['layerSet']][[thisFactor]]))[,c(1:3)], target_GR=chrom_bins, getBinCountOnly=TRUE)   # need to convert to GRanges to map to chrom
# some duplication of regions in the binstats, not sure what's going on
# count is not really same as coverage. 

# ? create temp GRanges and use this to calc coverages
wave.GR <-  GRanges(thisChrom, waveList[[k]][['layerSet']][[thisFactor]])
mcols(wave.GR) <- data.frame(score=width(wave.GR))    # set the 'score to be the width, can now use that as a 

widthAve <- avgByBin(xpr=data.frame( value=width(wave.GR)), featureData= as.data.frame( wave.GR)[,c(1:3)], target_GR=chrom_bins) 


union(GRanges("chrX", IRanges(c(300, 400), c(350, 500)), mcols=data.frame(value1=c(1,2))),
      GRanges("chrX", IRanges(c(200, 400), c(250, 500)), mcols=data.frame(value2=c(1,2))))
# Nope!  Discards mcols

ov <- findOverlaps(chrom_bins, wave.GR)

# aggregate(wave.GR, by=chrom_bins, FUN=function(x) sum(width(x)))   # NOPE!

aggregate(wave.GR, by= queryHits(ov), FUN=sum)

#aggregate(waveList[[k]][['layerSet']][[thisFactor]], by= queryHits(ov), FUN=sum)


wave.IR <- IRanges(start(wave.GR) , end(wave.GR))
chromBins.iR <- IRanges(start(chrom_bins), end(chrom_bins))
aggregate(wave.IR, by=chromBins.iR, FUN=length) 

wave.start <- GRanges(thisChrom, IRanges(start=start(waveList[[k]][['layerSet']][[thisFactor]]), end=start(waveList[[k]][['layerSet']][[thisFactor]])), mcols=data.frame(score=width(waveList[[k]][['layerSet']][[thisFactor]])))

ov <- findOverlaps(chrom_bins, wave.start)   # use start so each can only overlap one bin

aggregate(mcols(wave.start)[,1], by= queryHits(ov), FUN=sum)
by(mcols(wave.start)[,1], queryHits(ov), FUN=sum)

indexCounts <- tapply(mcols(wave.start)[,1], queryHits(ov), FUN=sum)

chrom.data <- as.data.frame(chrom_bins)
for(k in c(2,4,8,16,32,64))  {
  wave.start <- GRanges(thisChrom, IRanges(start=start(waveList[[k]][['layerSet']][[thisFactor]]), end=start(waveList[[k]][['layerSet']][[thisFactor]])), mcols=data.frame(score=width(waveList[[k]][['layerSet']][[thisFactor]])))
  ov <- findOverlaps(chrom_bins, wave.start)   # use start so each can only overlap one bin
  indexCounts <- tapply(mcols(wave.start)[,1], queryHits(ov), FUN=sum)
  chrom.data[,paste0("wave.", k)] <- 0
  chrom.data[names(indexCounts),paste0("wave.", k)] <- indexCounts
}

#plot(chrom.data$start, chrom.data$wave.100)




# created multibin GR object with metacols for each iteration of waves (or each 10th).
# ALT take doublings: wave 2,4,8,16,32, 64 etc.

# can is use avgByBin() ?  possibly not.

# do i want sum of regions, sum of widths in bin or coverage()?


#plot multibin as lines on ideogram, repeat for different marks. 







 