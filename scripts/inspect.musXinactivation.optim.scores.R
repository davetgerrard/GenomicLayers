

# 
# 
# 
# #setwd("C:/Users/Dave/HalfStarted/predictFromSequence/")
# 
# target.dir <- "results/mus_X_inactivation/mXi.BatchOpt.2016.11.11.dpsf/"
# target.files <- dir(target.dir, pattern="x.inactivation.mm9.10.2000")
# 
# resTable <- data.frame()
# 
# for(thisFile in target.files) {
#   load(paste(target.dir, thisFile, sep="/"))
#   
#   
#   
#   thisRow <- data.frame(file=thisFile, finalRound= length(optim.scores.vec), max.score=max(optim.scores.vec),
#                         final.stateWidth = keep.XFS$bf.spreadRep.1$mods$H3K27me3$stateWidth,
#                         final.offset = keep.XFS$bf.spreadRep.1$mods$H3K27me3$offset)
#   resTable <- rbind(resTable, thisRow)
# 
# }    
# 
# 
# write.table(resTable, file=paste(target.dir, "parameterResultsTable.tab", sep="/"), quote=F, row.names=F, sep="\t")



setwd("C:/Users/Dave/HalfStarted/predictFromSequence/")
target.dir <- "results/mus_X_inactivation/mXi.BatchOpt.2016.11.11.dpsf/"
scores <- read.delim("results/mus_X_inactivation/mXi.BatchOpt.2016.11.11.dpsf/parameterResultsTable.tab")
nrow(scores)
# starting values were stateWidth=200,offset=350
# mutation functions were:-
# mut.bf.spreadRep <-assignListNodeByCharacter(mut.bf.spreadRep, function() round(runif(1,min=50,max=400)), charRef="mods,H3K27me3,stateWidth")
# mut.bf.spreadRep <-assignListNodeByCharacter(mut.bf.spreadRep, function() round(runif(1,min=0,max=1000)), charRef="mods,H3K27me3,offset")

# how many runs reached 980 rounds?
hist(scores$finalRound, breaks=20)
mean(scores$finalRound)
sum(scores$finalRound == 980)    #75%

# expect mean statewidth 225
# expect mean offset 500

hist(scores$final.stateWidth, breaks=15)
mean(scores$final.stateWidth)
sum(scores$final.stateWidth == 200)   # very few have original value

hist(scores$final.offset, breaks=15)
mean(scores$final.offset)
sum(scores$final.offset == 350)     # very few have original value

# cannot distinguish this from mutation pressure.


hist(scores$max.score, breaks=15)
mean(scores$max.score)
# quite a tight distribution around .553
# I wonder how this would compare to a 100 runs of the initial factor set.
#     BUT not a fair comparison as these values are potentially the max values from a random distribution already.
#   Fairer comparison would be to run the keep.XFS object x100 and get distribution of scores.
#     (that would also be a better way to run the optimisation, but heavy on compute).












