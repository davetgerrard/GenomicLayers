
# trying to load a large number of pre-computed (known) TF motifs as pwm.

library(MotifDb)
length (MotifDb)
sort (table (values (MotifDb)$dataSource), decreasing=TRUE)


motifs <- query (MotifDb,  'sox')
motifs <- query(query (MotifDb,  'sox'), 'hsapiens')   # search terms nested
#hsapiens
show(motifs)
#export(motifs)
subset(motifs)[1:5]
motifs[1]
is(motifs[1])
motifs[1]@elementMetadata
motifs@elementMetadata   # note $bindingSequence
motifs@listData

query (MotifDb,  'nanog')
query(query (MotifDb,  'NANOG'), 'hsapiens') 




# read in some MEME format motifs from a txt file (encode website).
library(PWMEnrich)
motif.file <- "http://compbio.mit.edu/encode-motifs/motifs.txt"
motifs <- readMotifs(motif.file, remove.acc = FALSE)  # does not import correctly - non-standard file format.  DOES contain NANOG (several variants)

system.file(package="PWMEnrich", dir="extdata", file="example.transfac")




library(MotIV)



library(Biostrings)
data(HNF4alpha)
HNF4alpha




# iRegulon motif set
# downloaded huge set of motifs in cluster buster format (pfms) http://gbiomed.kuleuven.be/english/research/50000622/lcb/data-resources
# stored locally at 
motif.dir <- "C:/Users/Dave/data/iRegulon/motifCollection/singletons/"
  library(PWMEnrich)

fileName <- "taipale-NRTTAATNATTAACN-HNF1A-full.cb"
motif.file <- paste(motif.dir , fileName, sep="/")
#motifs <- readMotifs(motif.file) 

readCB <- function(fileName, pwm=TRUE)  {
  require(Biostrings)
  # *.cb files are individual motifs 
  # header line beginning with >
  # one line per bp with counts of ACGT at each position
  
  lines = readLines(fileName)
  
  mId <- sub(">", "", lines[1])
  #lines[2:length(lines)]
  lineList <- strsplit(lines[2:length(lines)], "\t")
  motifMatrix <- matrix(as.integer(unlist(lineList)), nrow=4, byrow = FALSE) ;row.names(motifMatrix) <- c("A", "C", "G" , "T")    # in horizontal form
  # motifMatrix <- matrix(as.integer(lineList)), ncol=4, byrow = TRUE ) 
  cs.mm <- colSums(motifMatrix)
  #if(!sd(cs.mm) == 0)  {  # for some reason, this matrix has columns of uneven size.
  #  # add pseudocounts to some columns
  #  max.col <- max(cs.mm)
  #  n.add <- max.col - cs.mm 
    
  #}
  if(pwm) {
  return(list(ID= mId, n.sites = min(cs.mm), pwm= prop.table(motifMatrix,2)))    # Temp fix with uneven colSums. TODO:investigate why colSUms uneven.
  } else {
    return(list(ID= mId, n.sites = min(cs.mm), pwm= motifMatrix)) 
  }
  # need to return as a motif object, retaining any metadata (e.g. name).
  #pwm <- PWM(motifMatrix)   # PWM requires that all columns sum to the same total. 
  # also PWM is perhaps more precise 
  
}


readFB <- function(fileName, pwm=TRUE)  {
  # read motifs from large factor book file (online at time of writing at http://compbio.mit.edu/encode-motifs/motifs.txt)
  lines <- readLines(fileName)
  id.index <- grep(">", lines)
  id.lines <- sub(">", "", lines[id.index])
  mList <- list()
  for(i in 1:length(id.index))  {
    
    mId <- id.lines[i]
    m.startLine <- id.index[i] +1 
    if(i == length(id.index)) {  # it's the last one
        m.endLine <- length(lines)
      } else { 
        m.endLine <- id.index[i+1] - 1
      }
    lineList <- strsplit(lines[m.startLine:m.endLine], " ")
    fullMatrix <- matrix(unlist(lineList), nrow=5, byrow = FALSE)
    row.names(fullMatrix) <- c("deg", "A", "C", "G" , "T")   
    motifMatrix <- fullMatrix[2:5,] 
    mode(motifMatrix) <- "numeric"   # easiest way to get char matrix to numeric.
    #return(list(ID= mId, n.sites = NA, pwm= motifMatrix)) 
    mList[[i]] <- list(ID= mId, n.sites = NA, pwm= motifMatrix)  # might be better to use mId as ref [[]] but needs to be unique.
  }
  return(mList)
  
}



fileName <- "http://compbio.mit.edu/encode-motifs/motifs.txt"

my.pwm <- readFB(fileName)

# search all the names
id.vec <- unlist(lapply(my.pwm, FUN=function(x)  x$ID))
nanog.index <- grep("NANOG", id.vec, ignore.case = T)
require(seqLogo)
#par(mfrow=c(2,3))
#for( i in nanog.index) {
 # seqLogo(my.pwm[[i]]$pwm)     # cannot draw several seqLogo on same plot
#}   

# To get multip seqLogos on one chart 
# https://support.bioconductor.org/p/35240/

#modify the seqLogo function 
mySeqLogo = seqLogo::seqLogo 
bad = (sapply( body(mySeqLogo), "==", "grid.newpage()") | 
         sapply( body(mySeqLogo), "==", "par(ask = FALSE)")) 
body(mySeqLogo)[bad] = NULL 

#function for generation of pwms 
#norm = function(x) scale(x, center=FALSE, scale=colSums(x)) 

#define dimensions 
ncol <- 2 
nrow <- 3 

#create plot 
grid.newpage() 
j <- 1    # index of nanog.index
for(row.i in 1:nrow){ 
  for(col.i in 1:ncol){ 
    j.n <- nanog.index[j]
    vp <- viewport(x = (col.i-1)/ncol, y = 1-(row.i-1)/nrow, w = 1/ncol, h = 1/nrow, just = c("left", "top")) 
    pushViewport(vp)
    pwm = my.pwm[[j.n]]$pwm
    mySeqLogo(pwm) 
    #grid.text(sprintf("Row %d, Column %d", row, col), x=0.5, y=0.9, just="top")
    upViewport() 
    j <- j +1
  } 
}

# ~ 2/3 of the iRegulon dataset have non-matched column totals. 

my.pwm <- readCB(motif.file)
#colSums(mm)
require(seqLogo)
seqLogo(my.pwm$pwm)
seqLogo(my.pwm$pwm, ic.scale=F)

# functions in Biostrings
round(my.pwm$pwm, 2)
maxWeights(my.pwm$pwm)
maxScore(my.pwm$pwm)
reverseComplement(my.pwm$pwm)


# could try using PFMatrix or PWMatrix classes from TFBSTools   # PFMatrix: Depends R (>= 3.2.2)  !
fileName <- "taipale-NRTTAATNATTAACN-HNF1A-full.cb"
motif.file <- paste(motif.dir , fileName, sep="/")
my.pwm <- readCB(motif.file, pwm=FALSE)
colSums(my.pwm$pwm)  # different numbers in columns
require(TFBSTools)
new.pfm  <- PFMatrix(ID=my.pwm$ID,   profileMatrix=my.pwm$pwm)   # Depends R (>= 3.2.2)  !



# DEVELOPMENT -----------------------

# # some of the .cb file matrices do not add up the same totals. 
# countMatch <- logical()
# for(thisFile in list.files(motif.dir))  {
#   motif.file <- paste(motif.dir , thisFile, sep="/")
#   mm <- readCB(motif.file)
#   countMatch[thisFile] <- sd(colSums(mm)) == 0
# }
# table(countMatch)

