devtools::install_github("hadley/devtools")
library(devtools)

install_github(repo="davetgerrard/GenomicLayers",build_vignettes = TRUE)
install_github(repo="davetgerrard/GenomicLayers",build_vignettes = TRUE, force=TRUE)

install_git("git://github.com/davetgerrard/GenomicLayers.git")


require(GenomicLayers)  # does not seem to load dependencies.
library(GenomicLayers)
?runLayerBinding
vignette(package="GenomicLayers")
createBindingFactor.DNA_motif("test", patternString = "ACTGGGCTA")  # works now.
vignette("GenomicLayersIntroduction")

# process the Roxygen2 annotations
library(devtools)
setwd("../GenomicLayers")
document()



# create a new vignette
library(devtools)
#use_vignette("GenomicLayersIntroduction")  # only do once to set up vignette for first time
#use_vignette("MouseX-inacivation") # only do once to set up vignette for first time


require(Biostrings)
build_vignettes()  #  did not load the dependencies by itself.  need require(Biostrings) 
