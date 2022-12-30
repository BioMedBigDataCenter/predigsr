
#' @title seurat workflow started with raw seurat count object
#'
#' @description  preprocess of usual seurat workflow
#'
#' @details
#'      The Workflow includes:
#'      (1) Normalization[Seurat::NormalizeData;FindVariableFeatures;ScaleData;RunPCA]
#'      (2) Batch Effect Removal[harmony::RunHarmony]
#'      (3) Clustering--[Seurat::RunUMAP;FindNeighbors;FindClusters;RunTSNE]
#' @param RawCountSeurat seurat object without "data" slot in the RNA assays
#'
#' @return
#' @import Seurat harmony
#' @export
#' @tests
#' alldata@meta.data$stim <- alldata@meta.data$donor
#' alldata <- get_seurat_workflow(alldata)
#' expect_equal(colnames(alldata@meta.data),colnames(newdata@meta.data))
get_seurat_workflow <- function(RawCountSeurat){
  alldata <- RawCountSeurat
  alldata <- alldata %>%
    Seurat::NormalizeData() %>%
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData(verbose = FALSE) %>%
    Seurat::RunPCA(pc.gene = alldata@var.genes,npcs = 20,verbose = FALSE)

  alldata <- alldata %>%
    harmony::RunHarmony("stim",plot_convergence = FALSE)
  alldata <- alldata %>%
    Seurat::RunUMAP(reduction="harmony",dim = 1:10) %>%
    Seurat::FindNeighbors(reduction="harmony",dims = 1:10) %>%
    Seurat::FindClusters(reduction="harmony",resolution = 1.5) %>%
    Seurat::RunTSNE()
  return(alldata)
}
