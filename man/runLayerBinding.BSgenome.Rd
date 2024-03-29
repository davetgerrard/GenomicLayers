% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runLayerBinding.BSgenome.R
\name{runLayerBinding.BSgenome}
\alias{runLayerBinding.BSgenome}
\title{Run the layer binding simulation on a BSgenome object}
\usage{
runLayerBinding.BSgenome(
  layerList,
  factorSet,
  iterations = 1,
  bf.abundances = rep(10000, length(factorSet)),
  watch.function = function(x) {
 },
  collect.stats = FALSE,
  keep.stats = TRUE,
  target.layer = 2,
  cache.layers = "LAYER.0",
  verbose = FALSE,
  ...
)
}
\arguments{
\item{layerList}{a \code{"Layerlist"} object containing a layerSet and other meta-data}

\item{factorSet}{a \code{"list"} of \code{"bindingFactor"} objects}

\item{iterations}{deprecated - do not use}

\item{bf.abundances}{10000  the quantities of each \code{"bindingFactor"} in \code{"factorSet"}}

\item{watch.function}{have this function execute during each iteration e.g. print something}

\item{collect.stats}{FALSE collect a table of stats each iteration}

\item{keep.stats}{TRUE  whether to retain data contained in $history.}

\item{target.layer}{NOT IMPLEMENTED}

\item{cache.layers}{store hits for these layers if they are immutable (e.g. LAYER.0 sequence) default= "LAYER.0". Set to NULL to prevent caching}

\item{verbose}{give more output}
}
\value{
\code{"LayerList"}
}
\description{
Simulate layer binding by the binding factor on the layerList
}
\examples{
# test to show whole genome layer binding 
# three factors are created, the third can only hit after the first has modified.

require(Biostrings)
require(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience


scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)

testFactor2 <- createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA")

testFactor3 <- createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                             mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

# check that a profile looking for 1 will not find any.  N.B this WILL bind AFTER testFactor2
testFactor4 <- createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
                                             mod.layers = c("LAYER.1", "LAYER.2"), mod.marks=c(0,1))

 
# now can match things genome wide. Need to run layerBinding and modification.

# need to have a factorSet, a list of bindingFactors

testFS <- list(testFactor2=testFactor2, testFactor3=testFactor3, testFactor4=testFactor4)

 
# with the above configuration, there are 41 possible sites across the genome, setting bf.abundances=30, restricts the number that are marked, so the number of potential sites reduces.
modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE, bf.abundances=30)

}
\seealso{
\code{\link{runLayerBinding}} \code{\link{matchBindingFactor.BSgenome}} \code{\link{modifyLayerByBindingFactor.BSgenome}} \code{\link{generateHitsCache}} \code{\link{createLayerSet.BSgenome}} \code{\link{createBindingFactor.DNA_motif}}
}
