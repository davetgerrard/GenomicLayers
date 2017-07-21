
# randomly select a value or its negative.
upDownFunc <- function(x)  {
  return(sample(c(x, -x), 1))
  
}

sapply(1:10, upDownFunc)

