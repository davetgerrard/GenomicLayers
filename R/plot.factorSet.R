
# show the profile and modification specificities for a set of binding factors.
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
  image(t(profile.layerMap), main="profile", axes=F, xlab="Layer") ; box()
  mtext(row.names(profile.layerMap), side=2, at=seq(0,1, length.out=nrow(profile.layerMap)), las=1)
  mtext(0:(ncol(profile.layerMap)-1), side=1, at=seq(0,1, length.out=ncol(profile.layerMap)), las=1)
  image(t(mod.layerMap), main="mods", axes=F, xlab="Layer"); box()
  mtext(row.names(mod.layerMap), side=2, at=seq(0,1, length.out=nrow(mod.layerMap)), las=1)
  mtext(1:ncol(mod.layerMap), side=1, at=seq(0,1, length.out=ncol(mod.layerMap)), las=1)
}
# plot.factorSet(factorSetRandom)
