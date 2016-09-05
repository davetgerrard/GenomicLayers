library(BSgenome.Hsapiens.UCSC.hg19)

genome <- BSgenome.Hsapiens.UCSC.hg19


tf.hits <- vmatchPattern("TTTCCCTAATC", genome, fixed=F)

tf.hits

gc()
