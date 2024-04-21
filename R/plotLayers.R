

plotLayers <- function(layerSet, layerNames, chrom, xlim) {
  
  
  # PLOTTING TO SHOW REAL TIME PROGRESS
  #xlim.gr <- c(3000000,3010000)  # mm10 chr19 shows unique banding for me1,me2, me3  - is an LTR
  #xlim.gr <- c(3762000,3772000) # kmt6b promoter and CpG island
  #xlim.gr <- c(47000000,47010000) 
  plot.new()
  plot.window(xlim, ylim= c(0,length(layerNames)))
  for(i in 1:length(layerNames)) {
    thisLayer <- layerNames[i]
    
    thisGR <- layerSet$layerSet[[thisLayer]]
    thisGR <- thisGR[seqnames(thisGR) == chrom]
    if(length(thisGR) > 0) {  
      rect(start(thisGR)-0.5, i-0.9 , end(thisGR)+0.5, i, col="black")
    }
  }
  
  axis(1)
  mtext(layerNames, side=2, at = 1:length(layerNames), las=2)
  # END OF PLOTTING
  
}
