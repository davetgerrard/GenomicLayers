
stopifnot(interactive())

gc()
system.time(
  results10 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland, verbose=T, max.window=10000000)
)

system.time(
  results20 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland, verbose=T, max.window=20000000)
)


system.time(
  results30 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland, verbose=T, max.window=30000000)
)

system.time(
  results40 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland, verbose=T, max.window=40000000)
)

system.time(
  results50 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland, verbose=T, max.window=50000000)
)

system.time(
  results60 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland, verbose=T, max.window=60000000)
)

system.time(
  results85 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland, verbose=T, max.window=85000000)   # fastest for 166M mus chrX
)

system.time(
  results100 <- matchBindingFactor(layerSet.X, bindingFactor = bf.CpGisland, verbose=T, max.window=100000000)
)