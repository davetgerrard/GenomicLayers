# mutateFactorSet
# description
# Parameters:-
# factorSet,
# layerSet ,
# mut_type="subRandomFactor",
# n.muts=1,
# verbose=FALSE,
# test.layer0.binding=FALSE,
# type_list = c("subRandomFactor" , "duplicate", "switch"),   "insert", "delete", "translocate", "move"
# prob=rep(1/length(type_list), length(type_list))
# fix.set.size= TRUE,   Keep the factor set the same size.
# name.prefix=""    give new factors a new name beginning with this name
#
# TODO re-write this so that a mixed set of mutation types can be generated.
mutateFactorSet <- function(factorSet, layerSet , mut_type="subRandomFactor", n.muts=1, verbose=FALSE, test.layer0.binding=FALSE,
                            type_list = c("subRandomFactor" , "duplicate", "switch"), prob=rep(1/length(type_list), length(type_list)), fix.set.size= TRUE, name.prefix=""){

  newFactorSet <- factorSet
  if(mut_type == "random")  {
   mut_type <- sample( type_list, 1, prob=prob)
  }


  if(mut_type== "duplicate")  {
    for(j in 1:n.muts)  {
      copy <- sample(1:length(newFactorSet), 1)  # which factor to duplicate
      index <- sample(1:length(newFactorSet), 1)  # what position to insert
      position  <- sample(c("before", "after"), 1)
      newName <- paste("dup", copy, sep=".")
      newFactor <- newFactorSet[[copy]]
      newFactor$name <-  newName
      # first join the newFactor
      append.index <- ifelse(position=="after", index, index-1)   # subtract one if "before"
      newFactorSet <- append(newFactorSet, list(newFactor), after=append.index)
      names(newFactorSet)[append.index+1] <- newName

      if(verbose)  {
        print(paste("Factor", copy,  "duplicated to ", position , "position", index ))
      }
    }

  }

  if(mut_type== "insert")  {
    for(j in 1:n.muts)  {
    index <- sample(1:length(newFactorSet), 1)  # what position to insert
    position  <- sample(c("before", "after"), 1)
    #for(i in index) {
    newName <- paste("ins", index, sep=".")
    newType <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"),1, replace=T)
    newFactor <- createRandomBindingFactor(newName,layerSet, type=newType, test.layer0.binding=test.layer0.binding, test.mismatch.rate=.1, verbose=verbose )
    #if(verbose) print(newFactor)
    # first join the newFactor
    append.index <- ifelse(position=="after", index, index-1)   # subtract one if "before"
    #fS.l <- length(factorSet)
    #append.vec <- append(1:fS.l, fS.l+1, after=append.index)
    newFactorSet <- append(newFactorSet, list(newFactor), after=append.index)
    names(newFactorSet)[append.index+1] <- newName
    #newFactorSet <- newFactorSet[c(append.vec)]  # reorder

    if(verbose)  {
        print(paste("Factor inserted", position , "position", index ))
    }
    #}
    }
  }

  if(mut_type== "delete")  {
    index <- sample(1:length(factorSet), n.muts)
    #rem.names <-   # when rework factor set, will need to remove the corresponding abundances (or only keep those we keep.)
    newFactorSet <- newFactorSet[-c(index)]
    if(verbose)  {
      print(paste("Factor(s)", paste(index, collapse=",") ,"deleted"))
    }

  }

  if(mut_type== "subRandomFactor")  {
    index <- sample(1:length(factorSet), n.muts)
    for(i in index) {
      thisName <- names(factorSet)[i]
      newName <- ifelse(name.prefix == "", thisName, paste(name.prefix, i, sep="."))
      newFactorSet[[thisName]] <-  createRandomBindingFactor(newName,layerSet, type=factorSet[[i]]$type, test.layer0.binding=test.layer0.binding, test.mismatch.rate=.1, verbose=verbose )
      if(verbose)  {
        print(paste("Factor", i ,"substituted"))
      }
    }

  }

  if(mut_type =="switch")  {  # switch a factor with another from the set
    index <- sample(1:length(factorSet), n.muts)
    for(i in index) {
      partner <- sample(setdiff(1:length(factorSet), i), 1)
      new.index <- 1:length(factorSet)
      new.index[c(i, partner)] <- c(partner, i)
      newFactorSet <- newFactorSet[new.index]
      if(verbose)  {
        print(paste("Factors", i ,"and", partner, "switched"))
      }
    }

  }


  return(newFactorSet)
}