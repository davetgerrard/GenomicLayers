
# may need to replace 'values' with a list.
buildMutationSkeleton <- function(charVec, values=NA,split=",")  {
  newModel <- list()
  if(length(values) == 1) values  <- rep(values, length(charVec))
  
  for(j in 1:length(charVec))  {
    charRef <- charVec[j]
    char.levels <- unlist(strsplit(charRef, split=split))
    #newModel[[char.vec[1]]] <- NULL
    for(i in 1:length(char.levels)) {
      if(!is.null(newModel[[char.levels[1:i]]])) {
        next;
      } else {
        if(i == length(char.levels)) {
          newModel[[char.levels[1:i]]] <-values[j]
        } else {
          newModel[[char.levels[1:i]]] <- list()
        }
      }
    }
    
  }
  return(newModel)
}