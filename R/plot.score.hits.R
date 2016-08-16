# plot multiple score measures. Green for good, red for bad.
plot.score.hits <- function(query, target, methods=c("acc", "tpr","tnr", "ppv", "npv", "fpr", "fnr"),
                            colours=structure(c("green", "green", "green", "green", "green", "red", "red"),
                                              names=c("acc", "tpr","tnr", "ppv", "npv", "fpr", "fnr")))  {
  scores <- numeric()
  for(i  in 1:length(methods)) {
    thisMethod <- methods[i]
    scores[thisMethod] <- score.hits(query=query, target = target, method = get(methods[i]))
  }
  #return(scores)
  barplot(scores, ylim=c(0,1), col=colours[methods])
}

# plot.score.hits(query=modLayerSet$layerSet$LAYER.5, target = tss.IR)