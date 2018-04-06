

# run polycomb over whole mouse genome. 
# track % of each genome covered with each iteration.
# run for limited number of iterations.

library(Biostrings)
library(GenomicLayers)
library(BSgenome.Mmusculus.UCSC.mm9)
# for plotting after sim.
library(ggplot2)
library(reshape2) # install.packages("reshape2")

#genome <- BSgenome.Mmusculus.UCSC.mm9 

mm9LayerSet <- createLayerSet.BSgenome(genome=BSgenome.Mmusculus.UCSC.mm9,n.layers = 3,  layer.names=c("CpG_island", "PRC", "H3K27me3"), verbose=TRUE)




bf.PRC <- createBindingFactor.layer_region("PRC",  patternLength=200, mismatch.rate=0, 
                                           profile.layers = "CpG_island", profile.marks = 1,  
                                           mod.layers = "PRC", mod.marks=1, stateWidth=500)

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

XFS <- list(PRC =bf.PRC, CpGisland= bf.CpGisland, spreadRep=bf.spreadRep)

#n.waves <- 100
n.factors <- 20000
n.runs <- 1000

timeList <- list()

finalLayer <- mm9LayerSet
i <- 1
statsTrace <- data.frame()
coverTrace <- data.frame()
thisTime <- system.time(
  while(i <= n.runs) {
    print(i)
    finalLayer <- runLayerBinding.BSgenome(layerList=finalLayer, factorSet=XFS, verbose=TRUE, iterations=n.factors, collect.stats = TRUE)
    thisRow <- cbind(i, finalLayer$history)
    statsTrace <- rbind(statsTrace , thisRow)
    covL <- lapply(coverage(finalLayer$layerSet$LAYER.4), sum)   # get sum of marked regions on specific chromosome
    thisRow <- cbind(i, as.data.frame(covL))
    coverTrace <- rbind(coverTrace, thisRow)
    if( i %% 10 == 0)  {   # export tables every 100th iteration
      write.table(statsTrace, file=paste0("simpleSpread.", n.runs, ".statsTrace.tab"), sep="\t", quote=F, row.names=F)
      write.table(coverTrace, file=paste0("simpleSpread.", n.runs, ".coverTrace.tab"), sep="\t", quote=F, row.names=F)
    }
    i <- i +1
  }
)

thisTime

#library(ggplot2)
#library(reshape2) # install.packages("reshape2")

#R = data.frame(Delta = c(1,2), UE = c(1,1), RE = c(3.8, 2.4))
meltCover = melt(coverTrace, id = "i")
pdf("mm9.polycomb.chromosomeCoverageOverTime.pdf")
ggplot(meltCover, aes(x = i, y = value, group = variable, colour = variable)) +
  geom_line()
# trajectories of each chromosome differ (given fixed amount of active substance)
#  (perhaps should scale for chromosome length).
dev.off()

