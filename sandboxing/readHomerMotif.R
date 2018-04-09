

fileName <- "../data/homerMotifs/motifs/oct4.motif"



readHomerMotif <- function(fileName, pwm=TRUE)  {
  require(Biostrings)
  # *.cb files are individual motifs 
  # header line beginning with >
  # one line per bp with counts of ACGT at each position
  
  lines = readLines(fileName)
  
  #mId <- sub(">", "", lines[1])
  #lines[2:length(lines)]
  #lineList <- strsplit(lines[2:length(lines)], "\t")
  lineList <- strsplit(lines, "\t")
  mId <- sub(">", "", lineList[[1]][2])
  homerConsensus <- sub(">", "", lineList[[1]][1])
  motifMatrix <- matrix(as.numeric(unlist(lineList[-1])), nrow=4, byrow = FALSE) ;row.names(motifMatrix) <- c("A", "C", "G" , "T")    # in horizontal form
  # motifMatrix <- matrix(as.integer(lineList)), ncol=4, byrow = TRUE ) 
  cs.mm <- colSums(motifMatrix)
  #if(!sd(cs.mm) == 0)  {  # for some reason, this matrix has columns of uneven size.
  #  # add pseudocounts to some columns
  #  max.col <- max(cs.mm)
  #  n.add <- max.col - cs.mm 
  
  #}
  if(pwm) {
    return(list(ID= mId, n.sites = min(cs.mm), pwm= prop.table(motifMatrix,2)))    # Temp fix with uneven colSums. TODO:investigate why colSUms uneven.
  } else {
    return(list(ID= mId, n.sites = min(cs.mm), pwm= motifMatrix)) 
  }
  # need to return as a motif object, retaining any metadata (e.g. name).
  #pwm <- PWM(motifMatrix)   # PWM requires that all columns sum to the same total. 
  # also PWM is perhaps more precise 
  
}