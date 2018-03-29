
# offset SacCer simulation to show spread of marks along chromosome after initial landing points.

library(Biostrings)
library(GenomicLayers)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# 
# pattern.core <- "GAAAC"
# pattern.long <- "TGAAACR"  # R = puRine (A,G)
# IUPAC_CODE_MAP

genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience
# genome
# seqnames(genome)
# summary(genome)
# organism(genome)
# provider(genome) 
# 

tf.hits <- vmatchPattern(pattern.long, genome, fixed=F) 


scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)


testFactor <- createRandomBindingFactor(name="test.1", layerSet=scLayerSet$layerSet, type="DNA_motif")


results <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor)



# need to prescribe a simple binding factor, made a simple function to create one.

testFactor2 <- createBindingFactor.DNA_motif("tf2", patternString="TGGGCTA")  # ACTGGGCTA does not hit all chromosomes

results <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor2)

bf.motif <- createBindingFactor.DNA_motif("tf3", patternString="TGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                             mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results3 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = bf.motif)




# check that a profile looking for 1 will not find any.  N.B this WILL bind AFTER testFactor2
testFactor4 <- createBindingFactor.DNA_motif("tf4", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(1,0), 
                                             mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results4 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor4)


# non-sequence bindingfactor with offset 
# want a binding factor that creates mods either up or downstream of binding site. New parameter 'offset.method' to accept function as argument
upDownFunc <- function(x)  {
  return(sample(c(x, -x), 1))
  
}
bf.spread4 <- createBindingFactor.layer_region("offset", patternLength = 5, profile.layers="LAYER.4", profile.marks=1, 
                                                mod.layers="LAYER.4", mod.marks = 1,  stateWidth = 7, offset = 150, offset.method=upDownFunc)

results5 <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=bf.spread4)
# now can match things genome wide. Need to run layerBinding and modification.

# need to have a factorSet, a list of bindingFactors

testFS <- list( bf.motif=bf.motif, bf.spread4=bf.spread4)


# also need for modifyLayerByBindingFactor.Views to work on BSgenome and hits

mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=testFactor3)


modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE)
# 2016-09-05 

# with the above configuration, there are 41 possible sites across the genome, setting iterations=30, restricts the number that are marked, so the number of potential sites reduces.
modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE, iterations=30, collect.stats = TRUE)

n.runs <- 5000

finalLayer <- scLayerSet
i <- 1
statsTrace <- data.frame()
coverTrace <- data.frame()
while(i <= n.runs) {
  finalLayer <- runLayerBinding.BSgenome(layerList=finalLayer, factorSet=testFS, verbose=TRUE, iterations=30, collect.stats = TRUE)
  thisRow <- cbind(i, finalLayer$history)
  statsTrace <- rbind(statsTrace , thisRow)
  covL <- lapply(coverage(finalLayer$layerSet$LAYER.4), sum)   # get sum of marked regions on specific chromosome
  thisRow <- cbind(i, as.data.frame(covL))
  coverTrace <- rbind(coverTrace, thisRow)
  if( i %% 100 == 0)  {   # export tables every 100th iteration
    write.table(statsTrace, file=paste0("simpleSpread.", n.runs, ".statsTrace.tab"), sep="\t", quote=F, row.names=F)
    write.table(coverTrace, file=paste0("simpleSpread.", n.runs, ".coverTrace.tab"), sep="\t", quote=F, row.names=F)
  }
  i <- i +1
}





#plot(statsTrace$i, statsTrace$Coverage.LAYER.4)
#plot(x=coverTrace$i, y=coverTrace[,-1], )

library(ggplot2)
library(reshape2) # install.packages("reshape2")

#R = data.frame(Delta = c(1,2), UE = c(1,1), RE = c(3.8, 2.4))
meltCover = melt(coverTrace, id = "i")
pdf("chromosomeCoverageOverTime.pdf")
ggplot(meltCover, aes(x = i, y = value, group = variable, colour = variable)) +
  geom_line()
# trajectories of each chromosome differ (given fixed amount of active substance)
#  (perhaps should scale for chromosome length).
dev.off()