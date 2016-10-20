createBindingFactor.DNA_regexp <- function(name,  type="DNA_regexp", patternString="N",patternLength=0, profile.layers=NULL,profile.marks=NULL,
                                           mod.layers=NULL,mod.marks=NULL, stateWidth=patternLength,
                                      test.layer0.binding=FALSE, test.mismatch.rate=.1 , max.pattern.tries=1000, min.DM.length=2, min.DR.length=10, verbose=FALSE) {
  
  # patternLength will be variable for regular expressions. Need separate parameter for modLength and may become a vector or list with different lengths for each layer.
  #patternLength <- nchar(patternString)
  #TODO sort out how to define patternLength or re-write other functions to accomodate variable patternLength
  profileList <- list(LAYER.0=list(pattern=patternString , mismatch.rate=0, length=patternLength))
  
  
  if(length(profile.layers) >0) {
  for(i in 1:length(profile.layers)) {
    thisLayer <- profile.layers[i]
    profileList[[thisLayer]] <- list(pattern=profile.marks[i], mismatch.rate=0.1, length=patternLength)
  }
  }
  modList <- list()
  if(length(mod.layers) >0) {
  for(i in 1:length(mod.layers)) {
  #for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
    thisLayer <- mod.layers[i]
    modState <- mod.marks[i]
    modList[[thisLayer]] <- list(state=modState, stateWidth=stateWidth, offset=0, align="centre")   # TODO make stateWidth independent of patternLength
  }
  }
  bindingFactor <- list(name=name, type=type,
                        profile=profileList,
                        mods=modList)
  
  return(bindingFactor)
  
}



createBindingFactor.DNA_regexp("test", patternString="ACTGGGCTA")
