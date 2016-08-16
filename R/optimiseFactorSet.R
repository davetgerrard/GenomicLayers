
# optimiseFactorSet
# description
# Parameters:-
# layerList,
# factorSet,
# testing.function,
# target.layer,
# target.vec,
# n.iter=10,
# target.score=1,
# mut.rate=0.1,
# modsPerCycle=100,
# test.layer0.binding=FALSE,
# method="fast",
# logFile="",
# logCycle=10, record scores every 'logCycle' iterations
# maxNoChange=n.iter,
# verbose=FALSE,
# use.parallel=FALSE,
# n.cores=1
# max.factors=50    # max number of factors allowed. (when reached, removes "insert" and "duplicate" mut_types
optimiseFactorSet <- function(layerList, factorSet, testing.function, target.layer, target.vec, n.iter=10, target.score=1, mut.rate=0.1, modsPerCycle=100,
                              test.layer0.binding=FALSE, method="fast", logFile="",logCycle=10, maxNoChange=n.iter, verbose=FALSE, use.parallel=FALSE, n.cores=1,
                              max.factors=50, mut_type_list = c("subRandomFactor" ,  "switch", "duplicate" , "insert", "delete"))  {

  if(logFile != "")  write.table(cbind("iter", "best.score"), row.names=F, col.names=F, sep="\t", quote=F,file=logFile)

  target.count <- length(target.vec)
  target.coverage <- sum(width(target.vec))
  chrom.size <- nchar(layerList$layerSet$LAYER.0)
  currentFactorSet <- factorSet
  scores.vector <- numeric()
  print("Calculating initial scores")
  if(method =="fast") {
    print("using fast algorithm")
    initialModLayer <- runLayerBinding(layerList=layerList, factorSet = factorSet, iterations=modsPerCycle, verbose=verbose)
  } else {
    print("using slow algorithm")
    initialModLayer <- runLayerBinding(layerList=layerList, factorSet = factorSet, iterations=modsPerCycle, verbose=verbose)
  }

  initialScore <- testing.function(initialModLayer, targetLayer=target.layer, target.vec=target.vec)
  if(is.na(initialScore )) {
    stop("initialScore was NA")

  }

  currentBestScore <- initialScore
  if(logFile != "") write.table(cbind(0, currentBestScore), row.names=F, col.names=F, quote=F,sep="\t",file=logFile, append=TRUE)
  iSinceLastImprovement <- 0  # use to break from following loop if very many rounds without improvement

  if(use.parallel) {
    # attempt at parallel runs
    require(parallel)
    mc <- getOption("cl.cores", n.cores)
    if(verbose) print(paste("Cores to use:" ,mc))

    cl <- makeCluster(mc)
    # TEMP need to export all functions to the nodes but don't know what they are.
    func.vec <- sapply(ls(.GlobalEnv), FUN =function(x) class(get(x)))
    func.names <- names(func.vec[func.vec == "function"])
    if(verbose) print(func.names)
    # TODO package code to avoid above
    clusterExport(cl, c( func.names, "layerList",
                         "modsPerCycle", "verbose", "target.layer","target.vec"), envir = environment())
    n.tries <- 10

  }


  for(i in 1:n.iter) {
    # write to log file if appropriate
    if((i %% logCycle ) == 0  & logFile != "") {
      write.table(cbind(i, currentBestScore), row.names=F, col.names=F, quote=F,sep="\t",file=logFile, append=TRUE)
      # could also save the current best factor set at this point.
      save(currentFactorSet, file=paste(dirname(logFile), "/", "currentFactorSet.",i,".Rdata",sep=""))
    }
    if(iSinceLastImprovement > maxNoChange)  {
      print(paste("No improvement in ", iSinceLastImprovement, " iterations, exiting!", sep=""))
      break ;
    }

    if(length(currentFactorSet) >= max.factors)  {   # stop the factor set growing larger than this
      use_mut_types <- setdiff(mut_type_list, c("insert", "duplicate"))
      if(verbose)  print(paste("Factor set limited to max.factors: ", max.factors))
    } else {
      use_mut_types <- mut_type_list

    }

    better <- FALSE
    if(use.parallel) {
      # create n.cores new factorSet by mutating the currentSet.
      factorSet.list <- replicate(mc, mutateFactorSet(currentFactorSet, layerList$layerSet,type_list = use_mut_types, mut_type= "random", n.muts=floor(length(currentFactorSet) * mut.rate), verbose=verbose, test.layer0.binding=T, name.prefix=paste("m",i,sep="")), simplify=FALSE)

      # list(...) added when cluster not picking up layerList
      #n.tries <- 10  # VERY WEIRDLY , the below command sometimes fails on first call but works on second (or later).  TODO: fix this!
      tries=0
      while(tries <= n.tries)  {
        worked <- try(mod.list <- parLapply(cl, factorSet.list, fun=function(x) runLayerBinding(layerList=layerList, factorSet = x, iterations=modsPerCycle, verbose=verbose)), TRUE)
        tries <- tries +1
        if(class(worked) != "try-error") {
          if(verbose) print(paste("parallel mod list tries: " , tries))
          break;
        }
      }
      if(verbose) print(paste(tries, "of", n.tries, "attempted"))
      #clusterExport(cl, "mod.list", envir = environment())
      par.scores <- parSapply(cl, mod.list, FUN=function(x) testing.function(layerList=x, targetLayer=target.layer, target.vec=target.vec))
      print(paste("x", mc, sep=""))


      best.index <- which.max(par.scores)
      newFactorSet <- factorSet.list[[best.index]]
      newModLayer <- mod.list[[best.index]]
      newScore <- par.scores[best.index]
    } else {
      newFactorSet <- mutateFactorSet(currentFactorSet, layerList$layerSet,n.muts=floor(length(currentFactorSet) * mut.rate),type_list = use_mut_types, mut_type= "random",  verbose=verbose, test.layer0.binding=test.layer0.binding, name.prefix=paste("m",i,sep=""))
      if(method =="fast") {
        #print("using fast algorithm")
        newModLayer <- runLayerBinding(layerList=layerList, factorSet = newFactorSet, iterations=modsPerCycle, verbose=verbose)
      } else {
        #print("using slow algorithm")
        newModLayer <- runLayerBinding(layerList=layerList, factorSet = newFactorSet, iterations=modsPerCycle, verbose=verbose)
      }
      newScore <- testing.function(newModLayer, targetLayer=target.layer, target.vec=target.vec)

    }

    n.regions <- length(newModLayer$layerSet[[target.layer]])
    coverage <- sum(width(newModLayer$layerSet[[target.layer]]))
    regionsWithHit <- sum(overlapsAny(newModLayer$layerSet[[target.layer]], target.vec))
    targetsWithHit <- sum(overlapsAny(target.vec, newModLayer$layerSet[[target.layer]]))


    scores.vector[i] <- newScore
    print(paste("Round", i, ". Facors:", length(newFactorSet),". Marks on target layer:", n.regions, ", Coverage:", coverage, ", Regions with a hit:", regionsWithHit,
                ", Targets Hit:", targetsWithHit, ", Chrom size:", chrom.size, ", Target count:", target.count, ", Target coverage:", target.coverage))
    print(paste("Round", i, ". OldScore", currentBestScore, "NewScore", newScore, "Better?"))
    if(is.na(newScore)) {
      print("Newscore = NA, skipping to next")
      next ;
    }
    if(newScore > currentBestScore) {
      better <- TRUE
      currentBestScore <- newScore
      currentFactorSet <- newFactorSet
      iSinceLastImprovement <- 0
      print("Yes!")
    } else {
      iSinceLastImprovement <- iSinceLastImprovement + 1
      print("No!")
    }
    #print(paste("Round", i, "oldScore", currentBestScore, "newScore", newScore, "Better?",better))

  }
  if(use.parallel)   stopCluster(cl)
  currentFactorSet$optimScores <- scores.vector
  return(currentFactorSet)
}
