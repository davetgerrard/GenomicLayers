

# Simulate 2-stage X-inactivation in mouse
#environment, packages and functions
# OPTIMISATION of offset of spreading feature. 


# Rscript.exe mus.x-inactivationModel.fromConfig.R --config CONFIG-FILE.R
.libPaths(c(.libPaths(), "C:/Users/Dave/Documents/R/win-library/3.3"))  # required to run from command line and pick up local package library
# Load CONFIG --------------------------------
library(getopt)
### this file has been setup to use a config script specifying the input and output files. 

spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'config'  , 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt <- getopt(spec);
if ( is.null(opt$config ) ) { opt$config <- 'scripts/CONFIG.mus.x-inactivationModel.Optimise.2016-11-04.R'}


source(opt$config)

if(!file.exists(OUTPUTDIR))  dir.create(OUTPUTDIR)

# TODO - CpG BF finding CpGs, but results not very life-like. Check paper, is it GC% instead?
# TODO - better/best method to assess results?



#source('scripts/pfs.functions.R')

# TODO add in some reports of what was loaded by the CONFIG file and which genome we're using.

# run several waves of XFS to show spreading.  
# need to collect data each time to plot. it. 


# runLayerBinding, score, mutate, runagain, score, if better, update.
new.XFS <- XFS
keep.XFS <- XFS
new.bf.spreadRep <- bf.spreadRep
optim.scores.vec <- numeric()


# run through once to get base score
waveList <- list()
waveList[[1]] <- current.X <- layerList.X
for(i in 2:n.waves) {
  print(paste("Running wave", i))
  current.X <- runLayerBinding(current.X, factorSet = new.XFS, iterations = n.iters, verbose=T)
  #waveList[[i]] <- current.X
}
#current.X <- runLayerBinding(layerList.X, factorSet = new.XFS, iterations = n.iters, verbose=T)
layerSubset <- restrict(current.X[['layerSet']][[thisFactor]], start = window.start, end=window.end)
metric.vec <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
current.score <- cor(meanScores[no.index], metric.vec[no.index] , use="complete.obs", method="spearman")
optim.scores.vec[1] <- current.score
best.score <- current.score

for(j in 2:n.optCycles) {
  print(paste("Optimisation round", j))
  
  bf.spreadRep.mutated <- mutateBindingFactor(new.bf.spreadRep, mut.bf.spreadRep, n.muts=1,verbose=TRUE)
  new.XFS <- list(bf.CpGisland =bf.CpGisland , bf.PRC.1=bf.PRC, bf.spreadRep.1=bf.spreadRep.mutated, bf.PRC.2=bf.PRC, bf.spreadRep.2=bf.spreadRep.mutated)
  
  current.X <- layerList.X  # reset to base value.
  for(i in 2:n.waves) {
    print(paste("Running wave", i))
    current.X <- runLayerBinding(current.X, factorSet = new.XFS, iterations = n.iters, verbose=T)
    #waveList[[i]] <- current.X
  }
  waveList[[j]] <- current.X
  
  layerSubset <- restrict(current.X[['layerSet']][[thisFactor]], start = window.start, end=window.end)
  metric.vec <- windowWidth(layerSubset, starts=bin.starts, window.size=window.length, limit=window.end)
  current.score <- cor(meanScores[no.index], metric.vec[no.index] , use="complete.obs", method="spearman")
  optim.scores.vec[j] <- current.score
  
  # test if better
  print(paste("Current best score: ", best.score))
  print(paste("New score: ", current.score))
  if(current.score > best.score ) {
    
    print("better score!")  
    new.bf.spreadRep <- bf.spreadRep.mutated
    keep.XFS <- new.XFS
    best.score <- current.score
  } else {
    print("no improvement") 
  }
  if(j %% saveEvery == 0)  {
    save(keep.XFS, optim.scores.vec, waveList, file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves, ".", n.iters,".", j,".waveList.Rdata"))
  }
  if(interactive())  {
    plot(1:length(optim.scores.vec), optim.scores.vec)
    abline(h=best.score)
  }
}

save(keep.XFS, optim.scores.vec, waveList, file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves, ".", n.iters,".waveList.Rdata"))






 