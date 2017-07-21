

# Simulate 2-stage X-inactivation in mouse
#environment, packages and functions

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
if ( is.null(opt$config ) ) { opt$config <- 'scripts/CONFIG.mus.x-inactivationModel.2016-10-27.R'}


source(opt$config)

if(!file.exists(OUTPUTDIR))  dir.create(OUTPUTDIR)

# TODO - CpG BF finding CpGs, but results not very life-like. Check paper, is it GC% instead?
# TODO - better/best method to assess results?



#source('scripts/pfs.functions.R')

# TODO add in some reports of what was loaded by the CONFIG file and which genome we're using.

# run several waves of XFS to show spreading.  
# need to collect data each time to plot. it. 


waveList <- list()
waveList[[1]] <- current.X <- layerList.X
for(i in 2:n.waves) {
  print(paste("Running wave", i))
  current.X <- runLayerBinding(current.X, factorSet = XFS, iterations = n.iters, verbose=T)
  waveList[[i]] <- current.X
}


save(XFS, waveList, file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", n.waves, ".", n.iters,".waveList.Rdata"))

# 
## run further waves - uncomment to run, takes a fair while. Saves every 'saveEvery' waves.
# n.waves <- 100
# load(paste0(OUTPUTDIR, "x.inactivation.", n.waves, ".waveList.Rdata"))
current.X <- waveList[[n.waves]]
next.wave <- n.waves + 1

for(i in next.wave:further.waves) {
  print(paste("Running wave", i))
  current.X <- runLayerBinding(current.X, factorSet = XFS, iterations = n.iters, verbose=T)
  waveList[[i]] <- current.X
  if(i %% saveEvery == 0)  {
    save(XFS, waveList, file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", i, ".", n.iters,".waveList.Rdata"))
  }
}

save(XFS, waveList, file=paste0(OUTPUTDIR, "x.inactivation.",GENOME,".", i, ".", n.iters,".waveList.Rdata"))


print(paste("Finished simulation" , RUN_NAME))







 