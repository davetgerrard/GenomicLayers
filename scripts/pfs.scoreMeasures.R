

# Some measures of specificity, sensitivity etc to use as scoring in optimisation.
###https://en.wikipedia.org/wiki/Sensitivity_and_specificity


acc <- function(TP, TN, FP, FN) { # accuracy
  
  return( (TP + TN) / (TP +FP +FN +TN))
  
}


# acc(50, 50, 3, 20)

# not suitable for genome data. Perhaps log it?
matthews.cor <- function(TP, TN, FP, FN) {
  return((TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) )
}

# matthews.cor(50, 50, 3, 20)


tpr <- function(TP, TN, FP, FN) {   #  true positive rate, TPR a.k.a.sensitivity
  return(TP/(TP+FN))
}
# tpr(50, 50, 3, 20)

tnr <- function(TP, TN, FP, FN) { #  true negative rate, TNR a.k.a. specificity
  return(TN/(TN+FP))
}
# tnr(50, 50, 3, 20)

ppv  <- function(TP, TN, FP, FN) { # positive predictive value, PPV a.k.a. precision
  return(TP/(TP+FP))
}
# ppv(50, 50, 3, 20)

npv <- function(TP, TN, FP, FN) { # negative predictive value, NPV
  return(TN/(TN+FN))
}
# npv(50, 50, 3, 20)


fpr <- function(TP, TN, FP, FN) { #false positive rate
  return(FP/(FP+TN))
}
# fpr(50, 50, 3, 20)

fnr <- function(TP, TN, FP, FN) { #false negative rate
  return(FN/(TP+FN))
}
# fnr(50, 50, 3, 20)


fdr <- function(TP, TN, FP, FN) { #false discovery rate
  return(FP/(TP+FP))
}

# fdr(50, 50, 3, 20)

# query and target are both IRanges objects
score.hits <- function(query, target, method=acc)  {
  FP <- sum(width(setdiff(query, target)))
  TP <- sum(width(intersect(query, target)))
  FN <- sum(width(setdiff(target, query)))
  TN <- sum(width(intersect(gaps(query), gaps(target))))
    
  score <- do.call(method, args=list(TP=TP, TN=TN, FP=FP, FN=FN))
  return(score)
}

# score.hits(query=modLayerSet$layerSet$LAYER.5, target=tss.IR)   # default method is accuracy
# score.hits(query=modLayerSet$layerSet$LAYER.5, target=tss.IR, method=fdr)
# score.hits(query=modLayerSet$layerSet$LAYER.5, target=tss.IR, method=matthews.cor)
# score.hits(query=modLayerSet$layerSet$LAYER.5, target=tss.IR, method=tnr)
# score.hits(query=modLayerSet$layerSet$LAYER.5, target=tss.IR, method=fnr)
