

assignListNodeByCharacter <-function(x, a, charRef=NULL, split=",")  {
  #charRef should be a single character string with successive names separated by the value of 'split'
  
  refVec <- unlist(strsplit(charRef, split=split))
  
  x[[refVec]] <- a
  return(x)
}