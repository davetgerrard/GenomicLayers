

#setwd("C:/Users/Dave/HalfStarted/predictFromSequence/")

target.dir <- "results/mus_X_inactivation/mXi.BatchOpt.2016.11.11.dpsf/"
target.files <- dir(target.dir, pattern="x.inactivation.mm9.10.2000")

resTable <- data.frame()

for(thisFile in target.files) {
  load(paste(target.dir, thisFile, sep="/"))
  
  
  
  thisRow <- data.frame(file=thisFile, finalRound= length(optim.scores.vec), max.score=max(optim.scores.vec),
                        final.stateWidth = keep.XFS$bf.spreadRep.1$mods$H3K27me3$stateWidth,
                        final.offset = keep.XFS$bf.spreadRep.1$mods$H3K27me3$offset)
  resTable <- rbind(resTable, thisRow)

}    


write.table(resTable, file=paste(target.dir, "parameterResultsTable.tab", sep="/"), quote=F, row.names=F, sep="\t")

