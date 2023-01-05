#' Formatted printing for a set of binding factors
#'
#' Formatted printing for a set of binding factors
#'
#' @param factorSet the list of binding factors.
#' 
#' @seealso \code{\link{plot.factorSet}} \code{\link{createBindingFactor.DNA_motif}}
#' 
#' @return NULL
#'
#' @examples
#' testFactor2 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA")
#' 
#' testFactor3 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
#'                                              mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))
#' 
#' # check that a profile looking for 1 will not find any.  N.B this WILL bind AFTER testFactor2
#' testFactor4 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
#'                                              mod.layers = c("LAYER.1", "LAYER.2"), mod.marks=c(0,1))
#' 
#' testFactor5 <- createBindingFactor.layer_region("test5", patternLength = 150)
#' # now can match things genome wide. Need to run layerBinding and modification.
#' 
#' # need to have a factorSet, a list of bindingFactors
#' 
#' testFS <- list(testFactor2=testFactor2, testFactor3=testFactor3, testFactor4=testFactor4, testFactor5=testFactor5)
#' 
#' print.bfSet(testFS)
#'
#' @export
print.bfSet <- function(factorSet) {
  print("A list of binding factors:-")
  bf.names <- names(factorSet)
  types <- lapply(factorSet, FUN=function(x) x$type)
  dna.patterns <- unlist(lapply(factorSet, FUN=function(x) as.character(x$profile$LAYER.0$pattern)))
  pattern.lengths <- unlist(lapply(factorSet, FUN=function(x) as.character(x$profile$LAYER.0$length)))
  layer.patterns <- unlist(lapply(lapply(factorSet, FUN=function(x) as.character(names(x$profile))), paste, collapse=","))
  layer.mods <- unlist(lapply(lapply(factorSet, FUN=function(x) as.character(names(x$mods))), paste, collapse=","))
  as.data.frame(cbind(bf.names, types,pattern.lengths, dna.patterns,  layer.patterns, layer.mods) )
}