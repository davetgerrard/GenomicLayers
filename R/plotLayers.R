#' Simple plotting function to monitor positions of features on a layerSet
#'
#' Simple plotting function to monitor positions of features on a layerSet
#'
#' @param layerSet the \code{"layerSet"} target
#' @param layerNames  which layers to plot
#' @param chrom    which chromosome to plot
#' @param xlim   a vector with start and end coordinates for the plotting window
#'
#' 
#'
#' @examples
#'\dontrun{
#' library(BSgenome.Scerevisiae.UCSC.sacCer3) 
#'
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#' 
#' sequences_to_keep <- setdiff(seqnames(genome), "chrM")  
#' 
#' genomeNuc <- keepBSgenomeSequences(genome, sequences_to_keep)
#' genomeNuc    # this should now still be a useable BSgenome object but with no mitochondrial chromosome.
#' # set up a LayerSet on the genome :  a list with link to genome and GRanges objects to store position information
#' # two  layers "H3K4me","H4K16ac" are specified here 
#' layerGenomeNuc <- createLayerSet.BSgenome(genome=genomeNuc, 
#'                                           layer.names=c( "H3K4me","H4K16ac"),
#'                                           n.layers=4)
#' # do something that creates data on the layers (e.g. run \code{"runLayerBinding"}                                          
#' layerGenomeNuc$layerSet[["H3K4me"]]    <- GRanges(seqnames="chrI", IRanges(start=c(400, 1000, 5000), end=c(600, 2000, 8000))                    
#'                                                   
#' plotLayers(layerGenomeNuc, layerNames=c( "H3K4me","H4K16ac"), chrom="chrI", xlim=c(1,10000)    
#'}
#' @import GenomicRanges
#' 
#' @export

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
