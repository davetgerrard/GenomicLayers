

getListNodeByCharacter <-function(x, charRef=NULL, split=",")  {
  
  nameParts <- unlist(strsplit(charRef, split=split))
  listObject <- x
  for(thisName in nameParts)  {   # should recurse down through list any length
    #print(thisName)
    listObject <- listObject[[thisName]]
    #print(x[[thisName]])
  }
  return(listObject)
}