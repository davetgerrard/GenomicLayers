#' Create a binding factor object to match a given pattern of layers, typically ignoring DNA
#'
#' Create a new binding factor based on a simple pattern of 
#' marks on others layers and \emph{may} (when used) set marks on other layers.
#'
#' @param name give the binding factor a name
#' @param type "layer_region"  to differentiate from other types
#' @param patternString  = NOT USED in this case
#' @param patternLength   [= 1] length of pattern
#' @param stateWidth the width of pattern to recognise on other layers, default is same as patternLength
#' @param profile.layers a vector of named layers to set as a match
#' @param profile.marks a vector of 0/1 to match the layers in profile.layers
#' @param mod.layers a vector of named layers to alter on a match
#' @param mod.marks a vector of 0/1 to set on the mod.layers
#' @param offset 0  integer value to indicate relative distance from pattern to apply modifications. Very simple. 
#' @param offset.method NULL   a \code{function} to apply to apply offset. MUST have parameter "n" that is used internally to represent the number of hits. 
#' @param offset.params NULL  a \code{list} of named parameters to pass to offset.method function
#' @param align "centre"
#' @param test.layer0.binding when creating, test if the DNA sequence has a match.
#' @param test.mismatch.rate proportion of mismatches to tolerate when testing
#' @param max.pattern.tries  NA
#' @param min.DM.length NA
#' @param min.DR.length NA
#' @param verbose set to TRUE for more output
#'
#' @return \code{"hits"}
#'
#' @examples
#' bf.LR <- createBindingFactor.layer_region("bf.LR", type="layer_region", 
#'          patternLength = 1, 
#'           stateWidth = 9, profile.layers = "LAYER.1",
#'           profile.marks = 0, mod.layers = "LAYER.1", mod.marks = 1)
#'          
#'  bf.LR1 <- createBindingFactor.layer_region("bf.LR1", type="layer_region",  patternLength = 1, patternString = "N",  profile.layers = "LAYER.1", profile.marks = 0, mod.layers = "LAYER.1", mod.marks = 1)
#'  bf.LR2 <- createBindingFactor.layer_region("bf.LR2", type="layer_region",  patternLength = 1, patternString = "N",  profile.layers = "LAYER.1", profile.marks = 0)  # profile but no mods
#'  bf.LR3 <- createBindingFactor.layer_region("bf.LR3", type="layer_region",  patternLength = 1, patternString = "N",   mod.layers = "LAYER.1", mod.marks = 1)  # no profile beyond LAYER.0 (genome)
#'  
#'  # example of using an offset.method to pass a function 
#'  upDownFuncRnorm <- function(n, offset.mean, offset.sd)  {
#'   y <- round(rnorm(n, mean=offset.mean, sd=offset.sd))
#'   z <- sample(c(1, -1), length(y), replace=T)  # random vector of 1,-1  to negate half the values
#'  
#'   return(round(y*z))
#'  }
#'  
#'  bf.LR4 <- createBindingFactor.layer_region("bf.LR4", type="layer_region",  
#'                     patternLength = 10,  profile.layers = "LAYER.1", profile.marks = 1, 
#'                     mod.layers = "LAYER.1", mod.marks = 1,
#'                    offset.method=upDownFuncRnorm, 
#'                    offset.params=list(offset.mean=50, offset.sd=15))
#'  
#' @export
createBindingFactor.layer_region <- function(name,  type="layer_region", patternLength=1, patternString=NULL,
                                             mismatch.rate=0, stateWidth=patternLength,
                                             profile.layers=NULL,  profile.marks=NULL,
                                             mod.layers=NULL,mod.marks=NULL,
                                             offset=0, offset.method=NULL, offset.params=NULL,
                                             align="centre",
                                            test.layer0.binding=FALSE, test.mismatch.rate=.1 , 
                                            max.pattern.tries=1000, min.DM.length=2, min.DR.length=10, verbose=FALSE) {
  
  # check input  
  stopifnot(exprs = {
    "profile.layers has non-unique names" = length(profile.layers) == length(unique(profile.layers))
    "mod.layers has non-unique names" = length(mod.layers) == length(unique(mod.layers))
  })  
  
  if(!is.null(offset.method)) stopifnot(exprs = { "offset.method if provided must be a function" = is.function(offset.method)   })
  
  #patternLength <- nchar(patternString)
  #profileList <- list(LAYER.0=list(pattern=DNAString(patternString) , mismatch.rate=0, length=patternLength))
  if(is.null(patternString)) {
    profileList <- list()
  } else {  # want to phase this out and have no DNA matching for this type
    profileList <- list(LAYER.0=list(pattern=DNAString(patternString) , mismatch.rate=mismatch.rate, length=patternLength))
  }
  
  if(length(profile.layers) >0) {  # there are layers to match beyond the sequence layer. Should always be true unless you want to match the whole chromosome.
    stopifnot("profile.marks does not match length of profile.layers" = length(profile.layers) == length(profile.marks))  
    for(i in 1:length(profile.layers)) {
      thisLayer <- profile.layers[i]
      profileList[[thisLayer]] <- list(pattern=profile.marks[i], mismatch.rate=mismatch.rate, length=patternLength)
    }
  }  
  
  modList <- list()
  if(length(mod.layers) >0) {
    stopifnot("mod.marks does not match length of mod.layers" = length(mod.layers) == length(mod.marks))
    for(i in 1:length(mod.layers)) {
      #for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
      thisLayer <- mod.layers[i]
      modState <- mod.marks[i]
      modList[[thisLayer]] <- list(state=modState, stateWidth=stateWidth, 
                                   offset=offset, offset.method=offset.method,
                                   offset.params=offset.params, align=align)   # TODO make stateWidth independent of patternLength
    }
  }
  
  bindingFactor <- list(name=name, type=type,
                        profile=profileList,
                        mods=modList)
  
  return(bindingFactor)
  
}



#bf.LR <- createBindingFactor.layer_region("layerBf", type="layer_region",  patternLength = 1, patternString = "N",  profile.layers = "LAYER.1", profile.marks = 0, mod.layers = "LAYER.1", mod.marks = 1)
#bf.LR <- createBindingFactor.layer_region("layerBf", type="layer_region",  patternLength = 1, patternString = "N",  profile.layers = "LAYER.1", profile.marks = 0, mod.layers = "LAYER.1", mod.marks = 1)
#bf.LR <- createBindingFactor.layer_region("layerBf", type="layer_region",  patternLength = 1, patternString = "N",  profile.layers = "LAYER.1", profile.marks = 0)  # profile but no mods
#bf.LR <- createBindingFactor.layer_region("layerBf", type="layer_region",  patternLength = 1, patternString = "N",   mod.layers = "LAYER.1", mod.marks = 1)  # no profile beyond LAYER.0 (genome)
