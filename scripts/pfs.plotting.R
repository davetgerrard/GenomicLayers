
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
colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs
print.bfSet(result[1:30])
load("data/HYDRA_runs/pfs_layer5_chr22_400bp_ppv/pfs_layer5_chr22_400bp_ppv.final.Rdata")
# plot.factorSet(factorSetRandom)
# plot.factorSet(result[1:30])
colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs
load("data/HYDRA_runs/pfs_layer5_chr22_400bp_tpr/pfs_layer5_chr22_400bp_tpr.final.Rdata")
# plot.factorSet(factorSetRandom)
# plot.factorSet(result[1:30])    # many more conversion of layer 5 to 1-state
colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs


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


# af.factorSet(result[1:(length(result)-1)])  # remember the last element is optimScores

system.time(modLayerSet <- runLayerBinding(layerList=layerList.5, factorSet = result[1:(length(result)-1)], verbose=TRUE, collect.stats = TRUE))  

modLayerSet$history

#raw.hits <- data.frame()
for(i in 1:(length(result)-1))  {
  
  matches <- matchBindingFactor(layerSet=layerList.5$layerSet, bindingFactor =result[[i]])
  thisRow <- data.frame(bf=names(result)[i], raw.hits=length(matches), raw.coverage=sum(width(matches)))
  if(i==1) {
    raw.hits <- thisRow
  } else {
    raw.hits <- rbind(raw.hits,thisRow)
  }
  
}

merge(raw.hits, modLayerSet$history)

