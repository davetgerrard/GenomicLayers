
# function(s) to describe the composition of a set of factors (a factorSet)

# perhaps could be used to describe the change in composition over an optimisation.




# plot.factorSet(factorSet)
# plot.factorSet(factorSetRandom)
# plot.factorSet(currentFactorSet)

load("data/HYDRA_runs/pfs_layer5_chr22_400bp_acc/pfs_layer5_chr22_400bp_acc.final.Rdata")
# plot.factorSet(factorSetRandom)
# plot.factorSet(result[1:30])    # many mods to silence (unmark) layer 5
colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs
print.bfSet(result[1:30])
load("data/HYDRA_runs/pfs_layer5_chr22_400bp_ppv/pfs_layer5_chr22_400bp_ppv.final.Rdata")
# plot.factorSet(factorSetRandom)
# plot.factorSet(result[1:30])
colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs
load("data/HYDRA_runs/pfs_layer5_chr22_400bp_tpr/pfs_layer5_chr22_400bp_tpr.final.Rdata")
# plot.factorSet(factorSetRandom)
# plot.factorSet(result[1:30])    # many more conversion of layer 5 to 1-state
colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs






# af.factorSet(result[1:(length(result)-1)])  # remember the last element is optimScores

system.time(modLayerSet <- runLayerBinding(layerList=layerList.5, factorSet = result[1:(length(result)-1)], verbose=TRUE, collect.stats = TRUE))  

modLayerSet$history

#raw.hits <- data.frame()
for(i in 1:(length(result)-1))  {
  
  matches <- matchBindingFactor(layerSet=layerList.5$layerSet, bindingFactor =result[[i]])
  thisRow <- data.frame(bf=names(result)[i], raw.hits=length(matches), raw.coverage=sum(width(matches)))
  if(i==1) {
    raw.hits <- thisRow
  } else {
    raw.hits <- rbind(raw.hits,thisRow)
  }
  
}

merge(raw.hits, modLayerSet$history)


load("data/HYDRA_runs/pfs_layer5_chr22_400bp_mutTest/pfs_layer5_chr22_400bp_mutTest.final.Rdata")
# plot.factorSet(factorSetRandom)
# n.factors <- length(result) -1
# plot.factorSet(result[1: n.factors])    # many more conversion of layer 5 to 1-state
colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs

load("data/HYDRA_runs/pfs_layer5_chr22_400bp_mutTest_200bf/pfs_layer5_chr22_400bp_mutTest_200bf.final.Rdata")
# plot.factorSet(factorSetRandom)
# n.factors <- length(result) -1
# plot.factorSet(result[1: n.factors])    # many more conversion of layer 5 to 1-state
colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs


# run against the full chromosome
library(BSgenome.Hsapiens.UCSC.hg19) # note the other packages being loaded.
genome <- BSgenome.Hsapiens.UCSC.hg19
thisChrom <- genome[["chr22"]] 

n.layers <- 5
layerSet.5 <- list(LAYER.0 = thisChrom)
for(i in 1:n.layers) {
  layerSet.5[[paste('LAYER.', i , sep="")]] <- IRanges()    # use IRanges to store state of layers. TODO limit to chrom length
}
layerList.5 <- list(layerSet=layerSet.5, history=NULL)
# run layerbinding with same number of mods as per optimisation
system.time(modLayerSet <- runLayerBinding(layerList=layerList.5, factorSet = result[1:(length(result)-1)], verbose=TRUE, collect.stats = TRUE, iterations=1000000,target.layer=5))  

modLayerSet$history

# plot the history showing the changing marked blocks and coverage on all the layers.
par(mfrow=c(1,2))
matplot(1:nrow(modLayerSet$history),modLayerSet$history[,c(4:8)], ylab="Coverage", xlab="Binding factor")
matplot(1:nrow(modLayerSet$history),modLayerSet$history[,c(9:13)], ylab="Number of blocks", xlab="Binding factor")
# DONE put into a function where the correct columns are found with grep



plot.layerBindingHistory(modLayerSet)

# what score? 
transcript.file <- "data/hg19.G19.chr22.transcript.gtf"
transcriptTable <- read.delim(transcript.file)
names(transcriptTable)[c(4,5,7)] <- c("txStart", "txEnd", "strand")
tss.positions <- ifelse(transcriptTable$strand == "+", transcriptTable$txStart, transcriptTable$txEnd)
tss.positions <- unique(tss.positions - base.0)    # only want unique values.
tss.positions <- na.omit(tss.positions)
tss.IR <- IRanges(start=tss.positions-199, end=tss.positions+200)   # this version of IRanges, resize not working right.


score.hits(modLayerSet$layerSet$LAYER.5, target = tss.IR, method = acc)






plot.score.hits(query=modLayerSet$layerSet$LAYER.5, target = tss.IR)




