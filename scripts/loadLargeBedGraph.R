




fileName <- "C:/Users/Dave/data/Simon_etal_2013_MusXist/bedGraphs/GSM1182888_d0.xist.mus.bedGraph.gz"


bg <- read.table(fileName, header=F, skip=1, sep= " ", nrows = 1000)


bg <- read.table(fileName, header=F, skip=1, sep= " ", nrows = 100000000)
