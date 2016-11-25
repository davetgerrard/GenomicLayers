
# attempts to create 'worse' versions of CpGisland BF

# calc score distribution (or just mean score?) for two new BF
# 1. like CpG island but without ANY sequence specificity (e.g. binds N)
# 2. like CpG island but with GC instead of CG.


# also need to get fast parallel scoring to work on cluster.


# 1. like CpG island but without ANY sequence specificity (e.g. binds N)
bf.CpGisland.N <- createBindingFactor.DNA_regexp("CpGisland", patternString="(CG.{0,20}){9}CG", patternLength=20,
                                               mod.layers = "CpG_island", mod.marks=1, stateWidth=200)


bf.CpGisland.N <- createBindingFactor.DNA_motif(name="CpGisland", type="DNA_motif", patternString="N", patternLength=108, 
                              mod.layers = "CpG_island", mod.marks=1, stateWidth=200)

#bf.CpGisland.N <- createBindingFactor.DNA_regexp("CpGisland", patternString="[NACGT]{108}", patternLength=20,
  #                                               mod.layers = "CpG_island", mod.marks=1, stateWidth=200)
#


results6 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland.N)

gc()

bf.GpCisland <- createBindingFactor.DNA_regexp("CpGisland", patternString="(GC.{0,20}){9}GC", patternLength=20,
                                                 mod.layers = "CpG_island", mod.marks=1, stateWidth=200)

results7 <- matchBindingFactor(layerSet.X, bindingFactor = bf.GpCisland)



# TODO check matchBindingFactor() to see if returning separate hits or reduced element.