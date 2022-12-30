#' reference_genes of digestive cells
#'
#' A data frame of genes used for cellular similarity corrleation calculation
#'
#' @format A data frame with 1 columns:
#' \describe{
#'   \item{Gene}{reference genes symbol or gene name}
#' }
#'
"reference_gene"

#' 0-1 value matrix of common cell types with their classical marker features
#'
#' A 0-1 value data.frame used for common cellular similarity corrleation calculation
#'
#' @format A data frame with 157 columns and 17 rows
#'
"data_matrix"

#' human-mouse gene conversion data frame
#'
#' A data frame with matched gene between human-hsa and mouse-mmu
#'
#' @format A data frame with 75483 columns and 4 rows
#'
#' @source \url{http://asia.ensembl.org/biomart/martview/65b916eaf1c6cf8f56e0ad08abf2a114}
#'
"trans"

#'
#'
#' A vector of 36 common color palette value used for ggplot2 visualization
#' @description Hexadecimal value
#' @source \url{https://www.jianshu.com/p/db0d0927e613}
"my36colors"
