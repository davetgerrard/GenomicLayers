
# THIS IS PROBABLY A DEAD END, FOR FIRST VIGNETTE, NEED SOMETHING MORE EPI-GENETIC...


## filter the SacCer gff file to TSS on chrI
# full gff from http://www.yeastgenome.org/download-data/sequence
# http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff

# gawk '$1 == "chrI" { print $0}' saccharomyces_cerevisiae.gff | awk '$3 == "gene" { print $0}' > saccharomyces_cerevisiae.chrI.genes.gff

?read.delim


ScChr1Genes <- read.table("c:/Temp/saccharomyces_cerevisiae.chrI.genes.gff", header=F)
head(ScChr1Genes)
nrow(ScChr1Genes)
promoterWidth <- 200
ScChr1Tss <- ifelse(ScChr1Genes$V7 == "+", ScChr1Genes$V4, ScChr1Genes$V5)
ScChr1Genes$PromStart <- ifelse(ScChr1Genes$V7 == "+", ScChr1Genes$V4 - promoterWidth,ScChr1Genes$V5)  
ScChr1Genes$PromEnd <- ifelse(ScChr1Genes$V7 == "+", ScChr1Genes$V4 ,ScChr1Genes$V5 + promoterWidth)

#ScChr1Genes.IR <- IRanges(ScChr1Genes$start, 

# http://genome.cshlp.org/content/25/7/1008.full
# "This results in transcription initiation between 40-120 bp downstream from the TATA element 
#   and typical core promoter lengths of 100-200 bp "

ScChr1Proms <- IRanges(start=ScChr1Genes$PromStart, end=ScChr1Genes$PromEnd)



modLayers$layerSet[['promoter']]

findOverlaps(modLayers$layerSet[['promoter']], ScChr1Proms)


# attempts to draw loci and predictinos. Not good.
library(utilsGerrardDT)

geneTable <- data.frame(chr=ScChr1Genes$V1, start=ScChr1Genes$V4, end=ScChr1Genes$V5, strand=ScChr1Genes$V7)
drawLoci(geneTable[25:50,])

TataHits <- scaleCoords(start(modLayers$layerSet[['promoter']]), input.range = c(1,max(geneTable$end)))
points(TataHits, y=rep(.5, length=length(TataHits)))

xlimits <- range(start(ScChr1Proms))  # show all
xlimits <- c(75000, 100000)
plot(x=start(ScChr1Proms),y=rep_len(c(.4, .5, .6), length(ScChr1Proms)), ylim=c(0,3), xlim=xlimits)
rect(xleft=start(ScChr1Proms), xright = end(ScChr1Proms), ybottom = 0, ytop = 3, col="grey70", border="grey70")
points(x=start(modLayers$layerSet[['promoter']]),y=rep_len(c(2.4, 2.5, 2.6), length(modLayers$layerSet[['promoter']])), pch=21,bg="black", col="black")
grid()

plot.new()
plot.window(xlim=xlimits, ylim=c(0,3))
rect(xleft=start(ScChr1Proms), xright = end(ScChr1Proms), ybottom = 0, ytop = 3, col="grey70", border="grey70")
points(x=start(modLayers$layerSet[['promoter']]),y=rep_len(c(2.4, 2.5, 2.6), length(modLayers$layerSet[['promoter']])), pch=21,bg="black", col="black")
randProms <- sample(1:max(start(ScChr1Proms)), length(modLayers$layerSet[['promoter']]))
points(x=randProms,y=rep_len(c(1.4, 1.5, 1.6), length(randProms)), pch=21,bg="red", col="red")

