

#' @title get_query_preprocess()
#' @description  Query Data Preprocessing Workflow for PreDigs
#'
#' @param seurat a query seurat object
#' @param species the species of query scRNAseq data ,i.e.,human or mouse
#' @param cellgroup the column name of cell type labels in the metadata data frame of each query seurat object
#' @param clusterlevel TRUE or FALSE, which TRUE means cell type annotation at cluster-level ,FALSE means annotation at single cell-level
#' @return
#' @export
#'
#' @examples
#' query_df <- get_query_preprocess(alldata,species = "human",cellgroup="cell_ontology_class",clusterlevel = TRUE)
#'
#' @tests
#' querydf <- get_query_preprocess(newdata,species = "human",cellgroup="celltype_match",clusterlevel = TRUE)
#' expect_equal(class(query_df),"data.frame")
get_query_preprocess <- function(seurat,species = "human",cellgroup,clusterlevel = TRUE){

  if (species == "human" & clusterlevel == TRUE){
    if (is.null(cellgroup)) {
      stop("Please provide the cellgroup label! i.e., cellgroup = seurat_clusters ")
    }
    print(cellgroup)
    # set cell type group in "celltype" column
    seurat@meta.data$celltype <- seurat@meta.data[,cellgroup]
    #print(colnames(seurat@meta.data))
    querydata_cor_df <- get_cor_pre_df(seurat)
  }

  if (species != "human" & clusterlevel == TRUE){
    if (is.null(cellgroup)) {
      stop("Please provide the cellgroup label! i.e., cellgroup = seurat_clusters")
    }
    # set cell type group in "celltype" column
    seurat@meta.data$celltype <- seurat@meta.data$cellgroup
    querydata_cor_df <- get_cor_pre_mmu_df(seurat)
  }

  if (species == "human" & clusterlevel == FALSE){
    querydata_cor_df <- get_cor_pre_sc_df(seurat)
  }

  if (species != "human" & clusterlevel == FALSE){
    querydata_cor_df <- get_cor_pre_sc_mmu_df(seurat)
  }
  ## the end
  return(querydata_cor_df)
}


#' @title get_cor_pre_df
#' @description  Get_cor_pre_df for human query data of healthy or cancer tissues at cluster-level
#'
#' @param seurat query data seurat object
#'
#' @return
#' @export
#'
get_cor_pre_df <- function(seurat){
  Anno_matrix_marker <- seurat@assays$RNA@data[which(rownames(seurat@assays$RNA@data) %in% reference_gene$Gene),]
  Anno_matrix_marker_df <- aggregate(as.data.frame(t(as.matrix(Anno_matrix_marker))),by=list(seurat@meta.data$celltype),FUN=mean)
  rownames(Anno_matrix_marker_df) <- Anno_matrix_marker_df$Group.1
  Anno_matrix_marker_df <- Anno_matrix_marker_df[,2:length(colnames(Anno_matrix_marker_df))]
  return(Anno_matrix_marker_df)
}

#' @title get_cor_pre_mmu_df
#' @description Get_cor_pre_mmu_df for mouse query data of healthy or cancer tissues at cluster-level
#'
#' @param seurat query data seurat object
#'
#' @return
#' @export
#'
get_cor_pre_mmu_df <- function(seurat){
  Anno_matrix_marker <- as.matrix(seurat@assays$RNA@data[,])
  Anno_matrix_marker <- Anno_matrix_marker[which(rownames(Anno_matrix_marker) %in% trans$Mouse.gene.name),]
  tmp_trans <- trans[which(trans$Mouse.gene.name %in% rownames(Anno_matrix_marker)),]
  tmp_trans <- tmp_trans[match(rownames(Anno_matrix_marker),tmp_trans$Mouse.gene.name),]
  Anno_matrix_marker_df <- aggregate.data.frame(Anno_matrix_marker,by=list(tmp_trans$Gene.name),median)
  rownames(Anno_matrix_marker_df) <- Anno_matrix_marker_df$Group.1
  Anno_matrix_marker_df <- Anno_matrix_marker_df[,2:length(colnames(Anno_matrix_marker_df))]
  Anno_matrix_marker_df <- aggregate(as.data.frame(t(as.matrix(Anno_matrix_marker_df))),by=list(seurat@meta.data$celltype),FUN=mean)
  rownames(Anno_matrix_marker_df) <- Anno_matrix_marker_df$Group.1
  Anno_matrix_marker_df <- Anno_matrix_marker_df[,2:length(colnames(Anno_matrix_marker_df))]
  return(Anno_matrix_marker_df)
}

#' @title get_cor_pre_sc_df
#' @description Get_cor_pre_sc_df for human query data of healthy or cancer tissues at single cell-level
#'
#' @param seurat query data seurat object
#'
#' @return
#' @export
#'
get_cor_pre_sc_df <- function(seurat){
  Anno_matrix_marker <- seurat@assays$RNA@data[which(rownames(seurat@assays$RNA@data) %in% reference_gene$Gene),]
  Anno_matrix_marker_df <- as.data.frame(t(as.matrix(Anno_matrix_marker)))
  return(Anno_matrix_marker_df)
}

#' @title get_cor_pre_sc_mmu_df
#' @description Get_cor_pre_sc_mmu_df for mouse query data of healthy or cancer tissues at single cell-level
#'
#' @param seurat query data seurat object
#'
#' @return
#' @export
#'
get_cor_pre_sc_mmu_df <- function(seurat){
  Anno_matrix_marker <- as.matrix(seurat@assays$RNA@data[,])
  Anno_matrix_marker <- Anno_matrix_marker[which(rownames(Anno_matrix_marker) %in% trans$Mouse.gene.name),]
  tmp_trans <- trans[which(trans$Mouse.gene.name %in% rownames(Anno_matrix_marker)),]
  tmp_trans <- tmp_trans[match(rownames(Anno_matrix_marker),tmp_trans$Mouse.gene.name),]
  Anno_matrix_marker_df <- aggregate.data.frame(Anno_matrix_marker,by=list(tmp_trans$Gene.name),median)
  rownames(Anno_matrix_marker_df) <- Anno_matrix_marker_df$Group.1
  Anno_matrix_marker_df <- Anno_matrix_marker_df[,2:length(colnames(Anno_matrix_marker_df))]
  Anno_matrix_marker_df <- as.data.frame(t(as.matrix(Anno_matrix_marker_df)))
  return(Anno_matrix_marker_df)
}
