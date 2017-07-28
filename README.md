## GenomicLayers - Simple, sequence-based simulation of epi-genomes.


Simulate changes to epigenetic states informed by genomic sequence.

Predict transcription (or any positional mark) from genome sequence by modelling the __factors__ that can bind to DNA and add or remove states. 

The raw DNA sequence cannot be altered but a series of __layers__ (binary vectors of same length of sequence) can have regions of any length switched (0 <-> 1). 

Some factors may also recognise patterns on the layers (e.g. regions in state 1) in addition to or instead of the underlying DNA sequence. 

The system models a series of very many sequential binding and modifying events. 

The ability of factors to bind changes through this series so that the number and positions of possible binding events for each factor change during the sequence. 

## Installation 

in R:-

	library(devtools)
	install_github(repo="davetgerrard/GenomicLayers",build_vignettes = TRUE)

The above may take several minutes and requires several dependencies. If it does not work, or you are in a hurry, leave out the 'build_vignettes'.

	install_github(repo="davetgerrard/GenomicLayers")


To view the introduction vignette, type

	vignette("GenomicLayersIntroduction")

To view what vignettes are available, type

	vignette(package="GenomicLayers")





