









pattern.core <- "GAAAC"
pattern.long <- "TGAAACR"  # R = puRine (A,G)
IUPAC_CODE_MAP

library(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience
genome
seqnames(genome)
summary(genome)
organism(genome)
provider(genome) 


tf.hits <- vmatchPattern(pattern.long, genome, fixed=F) 


scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)


testFactor <- createRandomBindingFactor(name="test.1", layerSet=scLayerSet$layerSet, type="DNA_motif")


results <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor)



# need to prescribe a simple binding factor, made a simple function to create one.

testFactor2 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA")

results <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor2)

# now can match things genome wide.







# will this run on human genome?  No, need much more memory to do human genome pattern matching.

library(BSgenome.Hsapiens.UCSC.hg19)

genome <- BSgenome.Hsapiens.UCSC.hg19

matchPattern("TTTCCCTAATC", genome, fixed=F)

tf.hits <- vmatchPattern("TTTCCCTAATC", genome, fixed=F) 
