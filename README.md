## GenomicLayers - Simple, sequence-based simulation of epi-genomes.


Simulate changes to epigenetic states informed by genomic sequence.

Predict transcription (or any positional mark) from genome sequence by modelling the __factors__ that can bind to DNA and add or remove states. 

The raw DNA sequence cannot be altered but a series of __layers__ (binary vectors of same length of sequence) can have regions of any length switched (0 <-> 1). 

Some factors may also recognise patterns on the layers (e.g. regions in state 1) in addition to or instead of the underlying DNA sequence. 

The system models a series of very many sequential binding and modifying events. 

The ability of factors to bind changes through this series so that the number and positions of possible binding events for each factor change during the sequence. 

## Installation 


GenomicLayers depends upon some [Bioconductor](https://www.bioconductor.org/) packages (BSGenome, GenomicRanges, Biostrings). One of the vignettes also depends on the BSgenome package BSgenome.Scerevisiae.UCSC.sacCer3 . If you install this first, it should install all the correct dependencies.

	if (!require("BiocManager", quietly = TRUE))
	  install.packages("BiocManager")
	BiocManager::install(version = "3.16")
	
	BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3") 

GenomicLayers can be installed direct from GitHub using the [devtools](https://github.com/hadley/devtools) package.  You will need the latest version of devtools.

Then to install GenomicLayers from Github:-

	library(devtools)
	install_github(repo="davetgerrard/GenomicLayers",build_vignettes = TRUE)

The above may take several minutes and requires several dependencies. If it does not work, or you are in a hurry, leave out the 'build_vignettes'.

	install_github(repo="davetgerrard/GenomicLayers")


To view the introduction vignette, type

	vignette("GenomicLayersIntroduction")

To view what vignettes are available, type

	vignette(package="GenomicLayers")





