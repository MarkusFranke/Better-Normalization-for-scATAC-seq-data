library(simATAC)

object <- newsimATACCount()
sim <- simATACSimulate(object, nCells = 100)
mtx <- simATACgetPeakByCell(sim, peak.num=5000)


object
sim