


mutateBindingFactor <- function(bf, mutSpectrum, n.muts=1, verbose=FALSE)  {
  new.bf <- bf
  if(verbose) print(paste(n.muts, "mutations requested"))
  #get vector of mutable list nodes
  vecList <- recurseListNames(mutSpectrum)
  
  if(verbose)  print(paste(length(vecList), "mutable characters"))
  
  for(i in 1:n.muts) {
    vec.index <- sample(1:length(vecList), 1)
    if(verbose)  print(paste("mutating", vecList[vec.index]))
    mut.call <- getListNodeByCharacter(mutSpectrum, vecList[vec.index])
    if(class(mut.call) == "function")  {
      #do.call(mut.call, args=list())
      new.bf <- assignListNodeByCharacter(new.bf, do.call(mut.call, args=list()), charRef=vecList[vec.index])  # currently not taking extra arguments
    } else {
      new.bf <- assignListNodeByCharacter(new.bf, mut.call, charRef=vecList[vec.index])
    }
  }
  return(new.bf)
}