
# attempts to create 'worse' versions of CpGisland BF

# calc score distribution (or just mean score?) for two new BF
# 1. like CpG island but without ANY sequence specificity (e.g. binds N)
# 2. like CpG island but with GC instead of CG.


# also need to get fast parallel scoring to work on cluster.








# runLayerBinding.parallel

# need to runLayerBinding in parallel to get distribution of scores.
# May be used on basic model on in optimisation.
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


# 1. like CpG island but without ANY sequence specificity (e.g. binds N)
bf.CpGisland.N <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){9}CG", patternLength=20,
                                               mod.layers = "CpG_island", mod.marks=1, stateWidth=200)


bf.CpGisland.N <- createBindingFactor.DNA_motif(name="CpGisland", type="DNA_motif", patternString="N", patternLength=108, 
                              mod.layers = "CpG_island", mod.marks=1, stateWidth=200)


# try setting patternString to empty string
bf.CpGisland.N <- createBindingFactor.DNA_motif(name="CpGisland", type="DNA_motif", patternString="", patternLength=108, profile.layers=NULL,
                                                mod.layers = "CpG_island", mod.marks=1, stateWidth=200)


#bf.CpGisland.N <- createBindingFactor.DNA_regexp("CpGisland", patternString="[NACGT]{108}", patternLength=20,
  #                                               mod.layers = "CpG_island", mod.marks=1, stateWidth=200)
#


results6 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland.N)

gc()

bf.GpCisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(GC.{0,20}){9}GC", patternLength=20,
                                                 mod.layers = "CpG_island", mod.marks=1, stateWidth=200)

results7 <- matchBindingFactor(layerSet.X, bindingFactor = bf.GpCisland)




bf.PRC <- createBindingFactor.layer_region("PRC",  patternLength=200, mismatch.rate=0, 
                                           profile.layers = "CpG_island", profile.marks = 1,  
                                           mod.layers = "PRC", mod.marks=1, stateWidth=500)

#bf.CpGisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){50}CG", patternLength=300, profile.layers="H3K27me3",profile.marks=0, mod.layers = "CpG_island", mod.marks=1)
# needed more seed areas, so took shorter CpG motifs.
# bf.CpGisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){9}CG", patternLength=20,
#                                                profile.layers="H3K27me3",profile.marks=0, 
#                                                mod.layers = "CpG_island", mod.marks=1, stateWidth=300)
#bf.CpGisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){9}CG", patternLength=20,
 #                                              mod.layers = "CpG_island", mod.marks=1, stateWidth=200)

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


n.waves <- 10
n.iters <- 2000

XFS <- list(bf.CpGisland =bf.CpGisland.N , bf.PRC.1=bf.PRC, bf.spreadRep.1=bf.spreadRep, bf.PRC.2=bf.PRC, bf.spreadRep.2=bf.spreadRep)

XFS.10 <- rep(XFS, times=n.waves)


# load comparative data for later testing
window.start <- 1
window.end <- length(thisChrom)

window.length <- 1*1000*1000

bin.starts <- seq(1, window.end, by=0.1*1000*1000)
# read in bedGraph Data from Simon et al.
fileName <- "C:/Users/Dave/data/Simon_etal_2013_MusXist/bedGraphs/GSM1182890_d3.xist.mus.bedGraph.gz"

bg <- read.table(fileName, header=F, skip=1, sep= " ")

nrow(bg)
names(bg) <- c("chrom", "start", "end", "value")
# only want chrX values.
bg.chrom <- subset(bg , chrom=="chrX")
nrow(bg.chrom)
bg.chrom.GR <- GRanges(bg.chrom$chrom, IRanges(start=bg.chrom$start, end=bg.chrom$end), value=bg.chrom$value)
bg.chrom.IR <- IRanges(start=bg.chrom$start, end=bg.chrom$end)
meanScores <- windowMean(bg.chrom.IR, y=mcols(bg.chrom.GR)$value, starts=bin.starts, window.size=window.length, limit=window.end)
no.index <- seq(30, length(meanScores), by=10) 
rm(bg, bg.chrom.IR, bg.chrom.GR)
gc()
thisFactor  <- "H3K27me3"
window.max <- 42*1000*1000
verbose <- TRUE


# combine runLayerBindng() and scoring to return a single value.

layerBindAndScore <- function(layerList, factorSet, iterations, thisFactor, window.start, window.end, window.length, 
                              bin.starts, meanScores, no.index )  {
  
  modLayer <- runLayerBinding( layerList=layerList, factorSet = factorSet, iterations =iterations)
  
  layerSubset <- restrict(modLayer[['layerSet']][[thisFactor]], start = window.start, end=window.end)
  metric.vec <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
  gc()
  return(cor(meanScores[no.index], metric.vec[no.index] , use="complete.obs", method="spearman"))
  
}


thisSCore <- layerBindAndScore(layerList=layerList.X, factorSet = XFS.10, iterations = n.iters*n.waves,
                               thisFactor=thisFactor, window.start=window.start, window.end=window.end,
                               window.length=window.length, bin.starts=bin.starts, meanScores=meanScores, 
                               no.index =no.index)  

# DONE check matchBindingFactor() to see if returning separate hits or reduced element.

require(parallel)
n.cores <- 4
mc <- getOption("cl.cores", n.cores)
if(verbose) print(paste("Cores to use:" ,mc))
cl <- makeCluster(mc)
# TEMP need to export all functions to the nodes but don't know what they are.
func.vec <- sapply(ls(.GlobalEnv), FUN =function(x) class(get(x)))
func.names <- names(func.vec[func.vec == "function"])
if(verbose) print(func.names)
# TODO package code to avoid above

n.waves <- 1
n.iters <- 2000 * 10   # to compensate for larger XFS.
clusterExport(cl, c( func.names, "layerList.X",
                     "XFS.10", "n.iters", "window.max", "n.waves",
                     "window.start", "window.end", "window.length",
                     "meanScores","no.index", "thisFactor", "bin.starts"), envir = environment())
n.tries <- 10     # this was because of a weird common failure on first try but not subsequent.


# for optimisation, will probably generate mutants and use parLapply (or parLapplyLB)  # don't export the factorset
# for distribution of results, will probably use clusterApply (or clusterApplyLB)   # export all functions and variables
system.time(
  resultsVec <- clusterApplyLB(cl, 1:30, fun=function(x) layerBindAndScore(layerList=layerList.X, factorSet = XFS.10, iterations = n.iters*n.waves,
                                                                          thisFactor=thisFactor, window.start=window.start, window.end=window.end,
                                                                          window.length=window.length, bin.starts=bin.starts, meanScores=meanScores, no.index =no.index))
)  # took XX on my PC.
#resultsSpread <- clusterCall(cl, 1:6, fun=function(x) runLayerBinding( layerList=layerList.X, factorSet = XFS, iterations = 10000))  #'works' but only gives one result per node.
stopCluster(cl)
gc()
resultsVec 

write(unlist(resultsVec), file="results/mus_X_inactivation/mXi.2016.12.09.null.CpG/scores.30.txt", ncolumns = 1)
png(file="results/mus_X_inactivation/mXi.2016.12.09.null.CpG/scores.30.bxp.png", res=150)
boxplot(unlist(resultsVec), ylim=c(-1,1))
dev.off()
