#' Formatted printing for a set of binding factors
#'
#' Formatted printing for a set of binding factors
#'
#' @param factorSet the list of binding factors.
#'
#' @return NULL
#'
#' @examples
#' x <- 1   # great!
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