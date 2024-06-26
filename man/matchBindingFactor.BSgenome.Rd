% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchBindingFactor.BSgenome.R
\name{matchBindingFactor.BSgenome}
\alias{matchBindingFactor.BSgenome}
\title{Find matches for a binding factor on a layer set containing a BSgenome sequence}
\usage{
matchBindingFactor.BSgenome(
  layerSet,
  bindingFactor,
  match.layers = names(bindingFactor$profile),
  clusterGap = 10,
  max.window = 1e+07,
  cache.layers = NULL,
  verbose = FALSE
)
}
\arguments{
\item{layerSet}{the \code{"layerSet"} target}

\item{bindingFactor}{the \code{"bindingFactor"} to match}

\item{match.layers}{restrict matches to only these named layers (default: all layers in names(bindingFactor$profile))}

\item{clusterGap}{=10 NOT IMPLEMENTED}

\item{max.window}{=10000000 on less powerful computers, break up the search into windows of this size.}

\item{cache.layers}{NULL which named layers to cache hits on (default NULL). Only use for fixed matches (e.g. to genome sequence "LAYER.0")}

\item{verbose}{output more information to the screen}
}
\value{
\code{"GRanges"}
}
\description{
Generate a list of matches for a binding factor against a layerSet object containing a BSgenome sequence. 
For this (.BSgenome) form, the \code{"Hits"} object returned will be a \code{"GRanges"} object
}
\examples{
require(Biostrings)
require(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience


scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)

testFactor3 <- createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                             mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

listOfHits <- matchBindingFactor.BSgenome(layerSet=scLayerSet, bindingFactor=testFactor3)


}
\seealso{
\code{\link{runLayerBinding}} \code{\link{runLayerBinding.BSgenome}} \code{\link{modifyLayerByBindingFactor.BSgenome}}
}
