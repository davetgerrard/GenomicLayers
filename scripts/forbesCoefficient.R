## Forbes coefficient
#https://academic.oup.com/bib/article/21/5/1523/5586919#230786692

## (N * TP ) / ((TP+FN)*(TP+FP))


forbes <- function(TP, FP, TN, FN, N)  {
  
  return((N * TP ) / ((TP+FN)*(TP+FP)))
}


#       pre
#  exp   T    F
#    T 450   30
#    F  50 1500

#+ve association, forbes = 3.8
contTable <- matrix(c(450, 30, 50, 1500), nrow=2, byrow=T, dimnames=list(exp=c("T", "F"), pre=c("T", "F")))
forbes(TP=contTable[1,1], FP=contTable[2,1], TN=contTable[2,2], FN=contTable[1,2], N=sum(contTable))

#perfect +ve association, forbes = 4.3   (1 + 1500/450)
contTable <- matrix(c(450, 0, 0, 1500), nrow=2, byrow=T, dimnames=list(exp=c("T", "F"), pre=c("T", "F")))
forbes(TP=contTable[1,1], FP=contTable[2,1], TN=contTable[2,2], FN=contTable[1,2], N=sum(contTable))

#perfect +ve association but bigger, forbes = 4.3   
contTable <- matrix(c(45000, 0, 0, 150000), nrow=2, byrow=T, dimnames=list(exp=c("T", "F"), pre=c("T", "F")))
forbes(TP=contTable[1,1], FP=contTable[2,1], TN=contTable[2,2], FN=contTable[1,2], N=sum(contTable))

#perfect  +ve association but very skewed, forbes = 3334.333  
contTable <- matrix(c(45, 0, 0, 150000), nrow=2, byrow=T, dimnames=list(exp=c("T", "F"), pre=c("T", "F")))
forbes(TP=contTable[1,1], FP=contTable[2,1], TN=contTable[2,2], FN=contTable[1,2], N=sum(contTable))

#no association, 50:50, forbes=1
contTable <- matrix(c(450, 450, 450, 450), nrow=2, byrow=T, dimnames=list(exp=c("T", "F"), pre=c("T", "F")))
forbes(TP=contTable[1,1], FP=contTable[2,1], TN=contTable[2,2], FN=contTable[1,2], N=sum(contTable))

#no association forbes=1
contTable <- matrix(c(10, 10, 450, 450), nrow=2, byrow=T, dimnames=list(exp=c("T", "F"), pre=c("T", "F")))
forbes(TP=contTable[1,1], FP=contTable[2,1], TN=contTable[2,2], FN=contTable[1,2], N=sum(contTable))

#neg association forbes=0.25
contTable <- matrix(c(10, 100, 200, 300), nrow=2, byrow=T, dimnames=list(exp=c("T", "F"), pre=c("T", "F")))
forbes(TP=contTable[1,1], FP=contTable[2,1], TN=contTable[2,2], FN=contTable[1,2], N=sum(contTable))


