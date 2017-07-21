windowWidth <- function(x, starts=NULL, window.size=NULL,  n.windows=length(starts), limit=max(starts)+1) {
  require(IRanges)
  # return the width of all features in in consecutive windows.
  # x should be an IRanges object
  
  
  
  width.vec <- integer()
  
  for(i in 1:n.windows)  {
    thisWindow.IR <- restrict(IRanges(start=starts[i], width=window.size), start = 1, end=limit)
    width.vec[i] <- sum(width(intersect(x, thisWindow.IR)))
    
  }
  return(width.vec)
  
}