createBindingFactor.DNA_motif <- function(name,  type="DNA_motif", patternString="N",
                                          patternLength = nchar(patternString), stateWidth=patternLength,
                                          profile.layers="LAYER.1",profile.marks=0,
                                          mod.layers="LAYER.1",mod.marks=1,
                                      test.layer0.binding=FALSE, test.mismatch.rate=.1 , max.pattern.tries=1000, min.DM.length=2, min.DR.length=10, verbose=FALSE) {
  
  
  profileList <- list(LAYER.0=list(pattern=DNAString(patternString) , mismatch.rate=0, length=patternLength))
  
  
  if(length(profile.layers) >0) {
  for(i in 1:length(profile.layers)) {
    thisLayer <- profile.layers[i]
    profileList[[thisLayer]] <- list(pattern=profile.marks[i], mismatch.rate=0.1, length=patternLength)
  }
  }
  modList <- list()
  for(i in 1:length(mod.layers)) {
  #for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
    thisLayer <- mod.layers[i]
    modState <- mod.marks[i]
    modList[[thisLayer]] <- list(state=modState, stateWidth=stateWidth, offset=0, align="centre")   #  make stateWidth independent of patternLength
  }
  
  bindingFactor <- list(name=name, type=type,
                        profile=profileList,
                        mods=modList)
  
  return(bindingFactor)
  
}



createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA")
