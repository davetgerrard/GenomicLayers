#' Customised matching for regular expression binding factors to BSgenomes
#'
#' Generate a list of matches for a binding factor against a BSgenome object. 
#' With proper forward and reverse patterns, will return hits from both strands.
#'
#' @param genome a \code{"BSgenome"} object
#' @param forRegExp see \code{\link{createBindingFactor.DNA_regexp}} 
#' @param revRegExp see \code{\link{createBindingFactor.DNA_regexp}} 
#'
#' @return \code{"Hits"}
#' 
#' @seealso \code{\link{createBindingFactor.DNA_regexp}} 
#'
#'
#' @import GenomicRanges
#' 
#' @export


matchRegExp.BSgenome <- function(genome, forRegExp, revRegExp)  {
  
  if(is.null(seqnames(genome))) stop("seqnames not set for genome. Please set using seqinfo() or seqnames()") 
     
     bsParams <- new("BSParams", X=genome, FUN=gregexpr)  # set up params for using bsapply
     grepResultBS_F <- bsapply(bsParams, pattern=forRegExp)  # run gregexpr over each chromosome separately
     grepResultBS_R <- bsapply(bsParams, pattern=revRegExp)  # run gregexpr over each chromosome separately
     all.hits <- GRanges(seqinfo=seqinfo(genome))   # empty GRanges to collate the results from different chroms
     for(chromName in seqnames(genome)) {  # cycle through chroms and convert matches to GRanges
       grepResult <- grepResultBS_F[[chromName]]
       if(grepResult[[1]][1] == -1 ) {  # no grep hits, do nothing
         #win.hits <- IRanges() 
       } else {
         win.hits <- GRanges(chromName, IRanges(start= as.integer(grepResult[[1]]), width = attr(grepResult[[1]], which="match.length", exact=TRUE)),seqinfo=seqinfo(genome))
       }
       all.hits <- c(all.hits, win.hits)
       grepResult <- grepResultBS_R[[chromName]]
       if(grepResult[[1]][1] == -1 ) {  # no grep hits, do nothing
         #win.hits <- IRanges() 
       } else {
         win.hits <- GRanges(chromName, IRanges(start= as.integer(grepResult[[1]]), width = attr(grepResult[[1]], which="match.length", exact=TRUE)),seqinfo=seqinfo(genome))
       }
       all.hits <- c(all.hits, win.hits)
     }
     return(all.hits)
}


#hits <- matchRegExp.BSgenome(genome=genomeSub, forRegExp="GGTGT(.{0,3})GGTGT", revRegExp="ACACC(.{0,3})ACACC")
#getSeq(genomeSub, hits)
