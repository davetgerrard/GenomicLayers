

# pfs CONFIG file for Mus X inactivation model
# specify details of simulation.

# run thus:-
# Rscript.exe mus.x-inactivationModel.fromConfig.R --config CONFIG-FILE.R --job_id XXX

#setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Mmusculus.UCSC.mm9)  # to compare with data from Simon et al.
#require(Biostrings)   # included in mm10 package dependencies
source('scripts/pfs.functions.R')
RUN_NAME <- "mXi.BatchOpt.2016.11.14.dpsf/"
OUTPUTDIR <- paste0("results/mus_X_inactivation/", RUN_NAME)
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

mut.bf.spreadRep <-buildMutationSkeleton(c("mods,H3K27me3,stateWidth", "mods,H3K27me3,offset" ), values=c(20, 50))
mut.bf.spreadRep <-assignListNodeByCharacter(mut.bf.spreadRep, function() round(runif(1,min=50,max=400)), charRef="mods,H3K27me3,stateWidth")
mut.bf.spreadRep <-assignListNodeByCharacter(mut.bf.spreadRep, function() round(runif(1,min=0,max=1000)), charRef="mods,H3K27me3,offset")




XFS <- list(bf.CpGisland =bf.CpGisland , bf.PRC.1=bf.PRC, bf.spreadRep.1=bf.spreadRep, bf.PRC.2=bf.PRC, bf.spreadRep.2=bf.spreadRep)



#n.waves <- 10
n.waves <- 10
n.iters <- 2000
saveEvery <- 10
further.waves <- 200
n.optCycles <- 1000

# to be used in plotting
thisFactor <- "H3K27me3"

window.start <- 1
window.end <- length(thisChrom)
window.length <- 1*1000*1000

bin.starts <- seq(1, window.end, by=0.1*1000*1000)
wavesToShow <- c(2,5,10,50,10,200)



# load comparative data

# read in bedGraph Data from Simon et al.
fileName <- "data/Simon_etal_2013_MusXist/bedGraphs/GSM1182890_d3.xist.mus.bedGraph.gz"

bg <- read.table(fileName, header=F, skip=1, sep= " ")

nrow(bg)
names(bg) <- c("chrom", "start", "end", "value")
# only want chrX values.
bg.chrom <- subset(bg , chrom=="chrX")
nrow(bg.chrom)
bg.chrom.GR <- GRanges(bg.chrom$chrom, IRanges(start=bg.chrom$start, end=bg.chrom$end), value=bg.chrom$value)
bg.chrom.IR <- IRanges(start=bg.chrom$start, end=bg.chrom$end)
meanScores <- windowMean(bg.chrom.IR, y=mcols(bg.chrom.GR)$value, starts=bin.starts, window.size=window.length, limit=window.end)
no.index <- seq(30, length(meanScores), by=10)   # non-overlapping windows.

