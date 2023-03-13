
# method to subset a BSgenome for testing.   Speeds up searches considerably. 

library(GenomicLayers)
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10

# hack from https://support.bioconductor.org/p/83588/  to keep only some chromosomes within a BSgenome object.
keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

#sequences_to_keep <- paste0("chr", c(1:20, "X", "Y"))
sequences_to_keep <- "chr17"
genomeSub <- keepBSgenomeSequences(genome, sequences_to_keep)
genomeSub

CGI<- createBindingFactor.DNA_regexp("CGI", patternString="(CG.{0,4}){3}CG", patternLength=20,
                                     mod.layers = "CpG_island", mod.marks=1, stateWidth=20)



#layerFullGenome <- createLayerSet.BSgenome(genome=genome, 
 #                                          layer.names=c("CpG_island",  "H3K27me1",  "H3K27me2","H3K27me3"),
#                                           n.layers=4)
layerPartialGenome <- createLayerSet.BSgenome(genome=genomeSub, 
                                              layer.names=c("CpG_island",  "H3K27me1",  "H3K27me2","H3K27me3"),
                                              n.layers=4)



gc40 <- createBindingFactor.DNA_regexp(name="gc40",
                                       patternString="")


bfSet <- list(CGI=CGI)   # binding factors must be in a list and named. Easiest to use each BFs name.
                                          

system.time(
  newLayerPartialGenome <- runLayerBinding.BSgenome(layerList=layerPartialGenome, factorSet=bfSet, iterations=100000) 
)


vmatchPattern("GCGCGCGGCGCG", subject=genomeSub,fixed=T)
vmatchPattern("GCGCGCGGCGCG", subject=genomeSub,fixed=T, max.mismatch = 2)
# to use the IUPAC codes (e.g. S for G/C), specify fixed =F. this then matches all the Ns at the telomeres.
vmatchPattern("GCGCGCGGCGCG", subject=genomeSub,fixed=F)  # here, no "S" but fixed =F creates many more hits..
vmatchPattern("GCGCGCGGCGCG", subject=genomeSub,fixed="subject")

vmatchPattern("GCGCGSGGCGCG", subject=genomeSub,fixed="subject")  # 13 regions 
vmatchPattern("GCGCGGGGCGCG", subject=genomeSub,fixed="subject")  # 10 hits G instead of S
vmatchPattern("GCGCGCGGCGCG", subject=genomeSub,fixed="subject")  # 3 hits C instead of S
vmatchPattern("GCGCGSGGCGCG", subject=genomeSub,fixed="pattern")  # matches many regions

vmatchPattern("SSSSSSSSSSSS", subject=genomeSub,fixed="subject")   # 51k regions
vmatchPattern( paste0(rep("S", 100), collapse=""), subject=genomeSub,fixed="subject")  # 100S has no hits
vmatchPattern( paste0(rep("S", 50), collapse=""), subject=genomeSub,fixed="subject")  #222 hits
vmatchPattern( paste0(rep("S", 100), collapse=""), subject=genomeSub,fixed="subject", max.mismatch = 10)  # 100S with 10 mismatches gives 1800 matches compared to 0 with
vmatchPattern( "SSSSSSSSSSSSSSSS", subject=genomeSub,fixed="subject", max.mismatch = 2)  # 100S with 10 mismatches gives 1800 matches compared to 0 with

# bf to match a given range of GC within 100bp
# finds regions of 100bp with 70-90bp of GC (S)
gc70_90 <- createBindingFactor.DNA_motif(name="gc70_90",
                                          patternString=paste0(rep("S", 100), collapse=""),
                                          max.mismatch= 30, 
                                          min.mismatch= 10,
                                          fixed="subject",
                                          profile.layers = NULL,
                                          profile.marks=NULL,
                                          mod.layers="H3K27me1", mod.marks = 1)

#this is identical but with fewer mismatches.
# mods to a different layer so that both can be run alongside each other.
gc80_100 <- createBindingFactor.DNA_motif(name="gc12",
                                      patternString=paste0(rep("S", 100), collapse=""),
                                      max.mismatch= 20, 
                                      min.mismatch= 0,
                                      fixed="subject",
                                      profile.layers = NULL,
                                      profile.marks=NULL,
                                      mod.layers="H3K27me2", mod.marks = 1)

bfSet <- list( gc70_90=gc70_90, gc80_100=gc80_100) 

newLayerPartialGenome <- runLayerBinding.BSgenome(layerList=layerPartialGenome, factorSet=bfSet, iterations=100000) 

length(newLayerPartialGenome$layerSet$H3K27me1)
length(newLayerPartialGenome$layerSet$H3K27me2)

#matchBindingFactor.BSgenome(layerSet=layerPartialGenome, bindingFactor=gc70_90)
#matchBindingFactor.BSgenome(layerSet=layerPartialGenome, bindingFactor=gc80_100)
