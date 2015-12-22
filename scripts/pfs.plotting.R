
# function(s) to describe the composition of a set of factors (a factorSet)

# perhaps could be used to describe the change in composition over an optimisation.


plot.factorSet <- function(factorSet)  {
  
  layerMap <- data.frame()
  
  profileLayers <- character()
  modLayers <- character()
  
  for(i in 1:length(factorSet)) {
    profileLayers <- sort(union(profileLayers, names(factorSet[[i]]$profile)))
    modLayers <-  sort(union(modLayers, names(factorSet[[i]]$mods)))
  }
  # could just set modLayers <- profileLayers ?
  
  profile.layerMap <- matrix(NA, ncol=length(profileLayers), nrow=length(factorSet), dimnames=list(factor=1:length(factorSet), layer=profileLayers ))
  mod.layerMap <- matrix(NA, ncol=length(modLayers), nrow=length(factorSet), dimnames=list(factor=1:length(factorSet), layer=modLayers ))
  
  for(i in 1:length(factorSet)) {
    
    
    
    for(thisLayer in names(factorSet[[i]]$profile)) {
      if(thisLayer == "LAYER.0") {
        profile.layerMap[i,thisLayer] <- 1
      } else {
        profile.layerMap[i,thisLayer] <- factorSet[[i]]$profile[[thisLayer]]$pattern
      }
      
    }
    
    
    for(thisLayer in names(factorSet[[i]]$mods)) {
      #mod.layerMap[i,thisLayer] <- TRUE  
      mod.layerMap[i,thisLayer] <- as.integer(factorSet[[i]]$mods[[thisLayer]]$state)
      
    }
    
  }
  
  par(mfrow=c(1,2))
  image(t(profile.layerMap), main="profile")
  image(t(mod.layerMap), main="mods")
  
}

# plot.factorSet(factorSet)
# plot.factorSet(factorSetRandom)
# plot.factorSet(currentFactorSet)

load("data/HYDRA_runs/pfs_layer5_chr22_400bp_acc/pfs_layer5_chr22_400bp_acc.final.Rdata")
# plot.factorSet(factorSetRandom)
# plot.factorSet(result[1:30])    # many mods to silence (unmark) layer 5
print.bfSet(result[1:30])
load("data/HYDRA_runs/pfs_layer5_chr22_400bp_ppv/pfs_layer5_chr22_400bp_ppv.final.Rdata")
# plot.factorSet(factorSetRandom)
# plot.factorSet(result[1:30])
load("data/HYDRA_runs/pfs_layer5_chr22_400bp_tpr/pfs_layer5_chr22_400bp_tpr.final.Rdata")
# plot.factorSet(factorSetRandom)
# plot.factorSet(result[1:30])    # many more conversion of layer 5 to 1-state


# create a matrix of nucleotide frequencies in the LAYER.0 patterns of the factorSet 
af.factorSet <- function(factorSet) {
  for(i in 1:length(factorSet)) {
    
    if(i == 1)  {
      af.table <-  alphabetFrequency(factorSet[[1]]$profile[['LAYER.0']]$pattern)
    }else {
      af.table <- rbind(af.table,alphabetFrequency( factorSet[[i]]$profile[['LAYER.0']]$pattern))
      
      
    }
    
    
  }
  return(af.table)
}


# af.factorSet(result)


