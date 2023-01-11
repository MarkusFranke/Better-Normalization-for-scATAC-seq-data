library(scran)
library(BiocParallel)

normScran <- function(data_mat, input_groups) {
  sizeFactors(
    computeSumFactors(
      SingleCellExperiment(
        list(counts = data_mat)),
      clusters = input_groups,
      min.mean = 0.1,
      BPPARAM = MulticoreParam()
    )
  )
}