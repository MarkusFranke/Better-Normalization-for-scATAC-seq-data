library(scran)
library(BiocParallel)
library(Seurat)
library(sctransform)

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

normSctransform <- function(adata) {
  seurat_obj <- as.Seurat(adata, counts = "X", data = NULL)
  seurat_obj <- RenameAssays(seurat_obj, originalexp = "RNA")
  res <- SCTransform(object = seurat_obj, method = "glmGamPoi", return.only.var.genes = FALSE, vst.flavor="v2")
  return(res)
}