# create a matrix of nucleotide frequencies in the LAYER.0 patterns of the factorSet
af.factorSet <- function(factorSet) {
  for(i in 1:length(factorSet)) {

    if(i == 1)  {
      af.table <-  alphabetFrequency(factorSet[[1]]$profile[['LAYER.0']]$pattern)
    }else {
      af.table <- rbind(af.table,alphabetFrequency( factorSet[[i]]$profile[['LAYER.0']]$pattern))


    }


  }
  return(af.table)
}
# colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs
