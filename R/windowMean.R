windowMean <- function(x, y,  starts=NULL, window.size=NULL,  n.windows=length(starts), 
                       window.func = mean, limit=max(starts)+1) {
  require(IRanges)
  # return the value of function on a certain score column in in consecutive windows.
  # x should be an IRanges object
  # y should be a vector of values of same lenght as x
  #     for several reasons, cannot yet use mcols() on x
  
  # there should be no self overlaps (otherwiise cannot make fair averages). 
  
  stopifnot( length(findOverlaps(x, drop.self=TRUE) ) ==0)   # should be no non-self overlaps.
  stopifnot(length(x) == length(y))
  
  mean.vec <- numeric()
  
  for(i in 1:n.windows)  {
    thisWindow.IR <- IRanges(start=starts[i], width=window.size)
    window.index <- overlapsAny(x, thisWindow.IR)
    
    # these two not currently used, was thinking to weight by number of regions or amount of space taken by regions within a window.
    # Will matter for uneven sized regions.
    window.coverage <- sum(width(x[window.index]))
    window.count <- sum(window.index)
    
    mean.vec[i] <- mean(y[window.index])
    
  }
  return(mean.vec)
  
}