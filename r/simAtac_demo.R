library(simATAC)

count <- getCountFromh5("r/GSE99172.snap")
dim(count)

object <- simATACEstimate(t(count))
a <- simATACget(object, "non.zero.pro")
setParameters(object, non.zero.pro=a)


object <- newsimATACCount()
sim <- simATACSimulate(object, nCells = 100)
mtx <- simATACgetPeakByCell(sim, peak.num=5000)


object
sim