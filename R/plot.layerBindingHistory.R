

plot.layerBindingHistory <- function(modLayerSet ) {

  coverage.cols <-  sort(grep("Coverage", names(modLayerSet$history), value=T))
  blockCount.cols <-  sort(grep("nBlocks", names(modLayerSet$history), value=T))
  par(mfrow=c(1,2))
  matplot(1:nrow(modLayerSet$history),modLayerSet$history[,coverage.cols], ylab="Coverage", xlab="Binding factor")
  matplot(1:nrow(modLayerSet$history),modLayerSet$history[,blockCount.cols], ylab="Number of blocks", xlab="Binding factor")

}

#plot.layerBindingHistory(modLayerSet)
