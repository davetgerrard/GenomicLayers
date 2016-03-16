# Load log of an optimisation run and plot some stats 

setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')



# Mutation test ------------------------------

runName <- "pfs_layer5_chr22_400bp_mutTest"
runLogFile <- "data/HYDRA_runs/pfs_layer5_chr22_400bp_mutTest/L5_c22.mutTest.400.o93515"
  
  
  logLines <- readLines(runLogFile)
  
  
  statLine.index <- grep("Coverage", logLines)
  scoreLine.index <- grep("OldScore", logLines)
  
  head(logLines[statLine.index])
  if(is.null(iters))  iters <- c(1 : length(statLine.index))
  
  # extract stats from each row
  
  #[1] "Round 545 . Facors: 50 . Marks on target layer: 433 , Coverage: 6900 , Regions with a hit: 63 , Targets Hit: 126 , Chrom size: 51304566 , Target count: 4070 , Target coverage: 1628000"
  #[1] "Round 545 . OldScore 0.970751174926647 NewScore 0.970690227249235 Better?"
  
  pattern <- ".*Facors: ([0-9]+) .*"
  factor.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))
  
  pattern <- ".*layer: ([0-9]+) .*"
  targetMarks.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))
  
  pattern <- ".*Coverage: ([0-9]+) .*"
  coverage.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))
  
  pattern <- ".*Regions with a hit: ([0-9]+) .*"
  regionsHit.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))
  
  pattern <- ".*Targets Hit: ([0-9]+) .*"
  targetsHit.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))
  
  pattern <- ".*OldScore ([0-9]\\.[0-9]+) .*"
  oldScores <- as.numeric(sub(pattern, "\\1", logLines[scoreLine.index]))
  
  pattern <- ".*NewScore ([0-9]\\.[0-9]+) .*"
  newScores <- as.numeric(sub(pattern, "\\1", logLines[scoreLine.index]))
  
  
  #plot(factor.count[1:100], type="l", ylim=c(0,80))
  #plot(factor.count, type="l", ylim=c(0,80))
  
  
  #plot(oldScores)
  
  #par(mfrow=c(2,1))
  #plot(oldScores)
  #plot(factor.count, type="l", ylim=c(0,80))
  
  png(paste("figures/", runName, ".statPlots.png", sep=""), height=1200, width=600, res=150)
  par(mfrow=c(6,1), mar=c(2,5,2,2))
  plot(oldScores[iters], ylab="Best score", xlab="", main=runName)
  plot(factor.count[iters], type="l", ylim=c(0,80), ylab="Number of factors", xlab="")
  plot(targetMarks.count[iters], ylab="Marks on target layer", xlab="")
  plot(coverage.count[iters], ylab="Coverage", xlab="")
  plot(regionsHit.count[iters], ylab="Regions with a hit", xlab="")
  plot(targetsHit.count[iters], ylab="Targets Hit", xlab="")
  dev.off()




# Accuracy test ----------------------
# note the late rise in coverage

runName <- "pfs_layer5_chr22_400bp_acc"
runLogFile <- "data/HYDRA_runs/pfs_layer5_chr22_400bp_acc/L5_c22.acc.400.o93111"
#png(paste("figures/", runName, "statPlots.png", sep=""), height=1200, width=600, res=150)
#plotLogs(runLogFile, iters=2:550)
#dev.off()

logLines <- readLines(runLogFile)


statLine.index <- grep("Coverage", logLines)
scoreLine.index <- grep("OldScore", logLines)

head(logLines[statLine.index])
 iters <- c(3 : length(statLine.index))   # for this run, skip the first iteration because it is so wildly different

# extract stats from each row

#[1] "Round 545 . Facors: 50 . Marks on target layer: 433 , Coverage: 6900 , Regions with a hit: 63 , Targets Hit: 126 , Chrom size: 51304566 , Target count: 4070 , Target coverage: 1628000"
#[1] "Round 545 . OldScore 0.970751174926647 NewScore 0.970690227249235 Better?"

pattern <- ".*Facors: ([0-9]+) .*"
factor.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))

pattern <- ".*layer: ([0-9]+) .*"
targetMarks.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))

pattern <- ".*Coverage: ([0-9]+) .*"
coverage.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))

pattern <- ".*Regions with a hit: ([0-9]+) .*"
regionsHit.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))

pattern <- ".*Targets Hit: ([0-9]+) .*"
targetsHit.count <- as.integer(sub(pattern, "\\1", logLines[statLine.index]))

pattern <- ".*OldScore ([0-9]\\.[0-9]+) .*"
oldScores <- as.numeric(sub(pattern, "\\1", logLines[scoreLine.index]))

pattern <- ".*NewScore ([0-9]\\.[0-9]+) .*"
newScores <- as.numeric(sub(pattern, "\\1", logLines[scoreLine.index]))


#plot(factor.count[1:100], type="l", ylim=c(0,80))
#plot(factor.count, type="l", ylim=c(0,80))


#plot(oldScores)

#par(mfrow=c(2,1))
#plot(oldScores)
#plot(factor.count, type="l", ylim=c(0,80))

png(paste("figures/", runName, ".statPlots.png", sep=""), height=1200, width=600, res=150)
par(mfrow=c(5,1), mar=c(2,5,2,2))
plot(oldScores[iters], ylab="Best score", xlab="", main=runName)
#plot(factor.count[iters], type="l", ylim=c(0,80), ylab="Number of factors", xlab="")
plot(targetMarks.count[iters], ylab="Marks on target layer", xlab="")
plot(coverage.count[iters], ylab="Coverage", xlab="")
plot(regionsHit.count[iters], ylab="Regions with a hit", xlab="")
plot(targetsHit.count[iters], ylab="Targets Hit", xlab="")
dev.off()



