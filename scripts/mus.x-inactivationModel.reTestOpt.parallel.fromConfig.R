

# Simulate 2-stage X-inactivation in mouse
#environment, packages and functions
# TESTING after optimisation. 
# Run as a SGE batch job (expecting numbered job-id on call).
# delete older results as we go.
# THIS VERSION: implements 16 core parallel duplicated runs at each stage, 

# Rscript.exe mus.x-inactivationModel.fromConfig.R --config CONFIG-FILE.R --job_id $SGE_TASK_ID
.libPaths(c(.libPaths(), "C:/Users/Dave/Documents/R/win-library/3.3"))  # required to run from command line and pick up local package library
# Load CONFIG --------------------------------
library(getopt)
### this file has been setup to use a config script specifying the input and output files. 

spec = matrix(c(
  'verbose', 'v', 0, "logical",
  'job_id', 'j', 1, "integer",
  'help'   , 'h', 0, "logical",
  'config'  , 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec);
if ( is.null(opt$config ) ) { opt$config <- 'C:/Users/Dave/HalfStarted/predictFromSequence/scripts/CONFIG.mus.x-inactivationModel.Optimise.2016-11-04.R'}
if (is.null(opt$verbose)) {opt$verbose <- FALSE }

job.id <- opt$job_id
verbose <- opt$verbose
source(opt$config)
print("make outputdir")
if(!file.exists(OUTPUTDIR))  dir.create(OUTPUTDIR)
par.outFile <- paste0(OUTPUTDIR,"/par.outFile.", opt$job_id, ".txt")  
# TODO - CpG BF finding CpGs, but results not very life-like. Check paper, is it GC% instead?
# TODO - better/best method to assess results?



#source('scripts/pfs.functions.R')

# TODO add in some reports of what was loaded by the CONFIG file and which genome we're using.

# run several waves of XFS to show spreading.  
# need to collect data each time to plot. it. 


# runLayerBinding, score, mutate, runagain, score, if better, update.
print("set  initial values")
n.cores <- 16
new.XFS <- XFS
keep.XFS <- XFS
XFS.10 <- rep(XFS, times=10)
new.bf.spreadRep <- bf.spreadRep
p.val.thres <- .10
#optim.scores.matrix <- numeric()

oldOutFile <- "notAFileName"    # to store the previous file name so that it can be deleted.

# run through once to get base score
#waveList <- list()
#waveList[[1]] <- current.X <- layerList.X
#for(i in 2:n.waves) {
#  print(paste("Running wave", i))
#  current.X <- runLayerBinding(current.X, factorSet = new.XFS, iterations = n.iters, verbose=T)
  #waveList[[i]] <- current.X
#}
#current.X <- runLayerBinding(layerList.X, factorSet = new.XFS, iterations = n.iters, verbose=T)
#layerSubset <- restrict(current.X[['layerSet']][[thisFactor]], start = window.start, end=window.end)
#metric.vec <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
#current.score <- cor(meanScores[no.index], metric.vec[no.index] , use="complete.obs", method="spearman")
#optim.scores.vec[1] <- current.score


### -> parallel

layerBindAndScore <- function(layerList, factorSet, iterations, thisFactor, window.start, window.end, window.length,
                              bin.starts, meanScores, no.index )  {
  print("begin layerBinding")
  modLayer <- runLayerBinding( layerList=layerList, factorSet = factorSet, iterations =iterations)
  print("calc score")
  layerSubset <- restrict(modLayer[['layerSet']][[thisFactor]], start = window.start, end=window.end)
  metric.vec <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
  gc()
  return(cor(meanScores[no.index], metric.vec[no.index] , use="complete.obs", method="spearman"))

}


print("Set up parallel cluster")

# now try on cluster
require(parallel)
mc <- getOption("cl.cores", n.cores)
print(paste("Cores available:", mc))
cl <- makeCluster(mc, outfile=par.outFile)
# TEMP need to export all functions to the nodes but don't know what they are.
func.vec <- sapply(ls(.GlobalEnv), FUN =function(x) class(get(x)))
func.names <- names(func.vec[func.vec == "function"])
if(verbose) print(func.names)
# TODO package code to avoid above
n.waves <- 1
n.iters <- 2000 * 10   # to compensate for larger XFS.
clusterExport(cl, c( func.names, "layerList.X",
                     "XFS.10", "n.iters", "window.max", "n.waves",
                     "window.start", "window.end", "window.length",
                     "meanScores","no.index", "thisFactor", "bin.starts"), envir = environment())
n.tries <- 10     # this was because of a weird common failure on first try but not subsequent.


# for optimisation, will probably generate mutants and use parLapply (or parLapplyLB)  # don't export the factorset
# for distribution of results, will probably use clusterApply (or clusterApplyLB)   # export all functions and variables
system.time(
  resultsVec <- clusterApplyLB(cl, 1:n.cores, fun=function(x) layerBindAndScore(layerList=layerList.X, factorSet = XFS.10, iterations = n.iters,
                                                                          thisFactor=thisFactor, window.start=window.start, window.end=window.end,
                                                                          window.length=window.length, bin.starts=bin.starts, meanScores=meanScores, no.index =no.index))
)  # took XX on my PC.
#resultsSpread <- clusterCall(cl, 1:6, fun=function(x) runLayerBinding( layerList=layerList.X, factorSet = XFS, iterations = 10000))  #'works' but only gives one result per node.
#stopCluster(cl)
gc()
#print(resultsVec)  # comes back as a list, 20 mins on my PC (4 runs on 3 cores)<


### <- parallel

#current.score <- min(unlist(resultsVec))
optim.scores <- best.scores <- current.scores <- matrix(sort(unlist(resultsVec)), ncol=n.cores)
print(best.scores)

#optim.scores.vec[1] <- current.score
#best.score <- current.scores



#for(j in 2:n.optCycles) {
  print("begin test")
  
  load(paste0(OUTPUTDIR, "x.inactivation.",GENOME,".",  n.iters,".",job.id, ".", n.optCycles,".scores.Rdata"))  # retrieve keep.XFS object
  #bf.spreadRep.mutated <- mutateBindingFactor(new.bf.spreadRep, mut.bf.spreadRep, n.muts=1,verbose=TRUE)
#  new.XFS <- list(bf.CpGisland =bf.CpGisland , bf.PRC.1=bf.PRC, bf.spreadRep.1=bf.spreadRep.mutated, bf.PRC.2=bf.PRC, bf.spreadRep.2=bf.spreadRep.mutated)
   XFS.10 <- rep(keep.XFS, times=10)

  #current.X <- layerList.X  # reset to base value.
  #for(i in 2:n.waves) {
  #  print(paste("Running wave", i))
  #  current.X <- runLayerBinding(current.X, factorSet = new.XFS, iterations = n.iters, verbose=T)
  #  #waveList[[i]] <- current.X
  #}
  #waveList[[j]] <- current.X
  
  #layerSubset <- restrict(current.X[['layerSet']][[thisFactor]], start = window.start, end=window.end)
  #metric.vec <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
  #current.score <- cor(meanScores[no.index], metric.vec[no.index] , use="complete.obs", method="spearman")
  
  clusterExport(cl, c( func.names, "layerList.X",
                     "XFS.10", "n.iters", "window.max", "n.waves",
                     "window.start", "window.end", "window.length",
                     "meanScores","no.index", "thisFactor", "bin.starts"), envir = environment())

   system.time(
  	resultsVec <- clusterApplyLB(cl, 1:n.cores, fun=function(x) layerBindAndScore(layerList=layerList.X, factorSet = XFS.10, iterations = n.iters,
                                                                          thisFactor=thisFactor, window.start=window.start, window.end=window.end,
                                                                          window.length=window.length, bin.starts=bin.starts, meanScores=meanScores, no.index =no.index))
   ) 
  #current.score <- min(unlist(resultsVec))
  # optim.scores.vec[j] <- current.score
 
  current.scores <- matrix(sort(unlist(resultsVec)), ncol=n.cores)
 
  optim.scores <- rbind(optim.scores, current.scores)

 
  # test if better
  cat("Current best scores: ", best.scores, "\n")
  cat("New scores: ", current.scores, "\n")
  wilcox.p <-  wilcox.test(current.scores,best.scores, alternative = "greater")$p.value
  print(paste("Wilcox p-value: ", wilcox.p))
   if( wilcox.p   <= p.val.thres ) {
    
    print("better scores!")  
    #new.bf.spreadRep <- bf.spreadRep.mutated
    #keep.XFS <- new.XFS
    #best.scores <- current.scores
  } else {
    print("no improvement") 
  }
  #if(j %% saveEvery == 0)  {
  #  outFile <- paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.iters,".",job.id, ".", j,".waveList.Rdata")
  #  save(keep.XFS, optim.scores, file=outFile)
  #  unlink(oldOutFile)   # remove the previous results file.
  #  oldOutFile <- outFile
  #}
  if(interactive())  {
    #plot(1:length(optim.scores.vec), optim.scores.vec)
    #abline(h=best.score)
  }


stopCluster(cl)
gc()


#save(keep.XFS, optim.scores.vec, waveList, file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves, ".", n.iters,".waveList.Rdata"))
outFile <- paste0(OUTPUTDIR, "x.inactivation.",GENOME,".",  n.iters,".",job.id, ".", j,".reTest.Rdata")
save(keep.XFS, optim.scores, wilcox.p, file=outFile)
#unlink(oldOutFile)   # remove the previous results file.
#oldOutFile <- outFile





 
