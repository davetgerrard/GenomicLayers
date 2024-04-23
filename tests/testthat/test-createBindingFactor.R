





#createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA" , profile.layers=c("bob", "bob")) 

test_that("function reject bad layer or mismatched layer names", {
  expect_error( createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA" , profile.layers=c("bob", "bob")) , "profile.layers has non-unique names")
  expect_error( createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA" , profile.layers=c("bob", "john")), "profile.marks does not match length of profile.layers") 
  expect_error( createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA" , profile.layers=c("bob", "john"), profile.marks=c(1)) ,  "profile.marks does not match length of profile.layers")
  expect_error( createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA" , mod.layers=c("bob", "bob")) , "mod.layers has non-unique names")
  expect_error( createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA" , mod.layers=c("bob", "john")), "mod.marks does not match length of mod.layers") 
  expect_error( createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA" , mod.layers=c("bob", "john"), mod.marks=c(1)) ,  "mod.marks does not match length of mod.layers")
})
