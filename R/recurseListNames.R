

# unlisted names are not suitable for parsing the list as different configurations could lead to the same unlisted name
# e.g.   a.b$c  and a$b.c both become a.b.c
recurseListNames <- function(x, stemVec=NULL, routeVec = NULL, sep=",")  {
  
  for(thisName in names(x)) {
    
    if(length(stemVec) < 1)  {  # first iteration
      newStem <- thisName
    }  else {
      newStem <- paste(stemVec, thisName, sep=sep)
    }
    #print(paste(stemVec, thisName))
    #stemList[[thisName]] <- c(stem, thisName)
    if( class(x[[thisName]]) == "list")  {
      routeVec <- recurseListNames(x[[thisName]], stemVec=newStem, routeVec=routeVec)
    } else {
      #stemVec <- paste(stemVec, thisName, sep=".") 
      #print(paste(stemVec, thisName))
      routeVec <- c(routeVec, newStem)
    }
    
  }
  return(routeVec)
}