

#' @title get_predigs_score_df()
#' @description  cell type prediction through predigs
#' @details calculate gene expression similarity matrix between query cells and reference cell types based on the multiple datasets integration
#' @param reference_cor_list a list of multi-reference datasets for cell type assignment from predigs
#' @param querydata_cor_list a list of multi-query datasets, and a query list with only one dataset is allowed
#' @param i the index value of querydata_cor_list, i.e., it could be 1,2,3,4,5 if the length of querydata_cor_list is equal to 5
#'
#' @importFrom reshape2 melt
#' @importFrom philentropy distance
#'
#' @return
#' @export
#'
#' @examples
#' data("ref_health_pancrea_cor_list")
#' querydata_list <- list(query_df)
#' query_predigs_score_list <- lapply(1:length(querydata_list),function(x){
#'     get_predigs_score_df(reference_cor_list = ref_health_pancrea_cor_list,
#'                                          querydata_cor_list = querydata_list,
#'                                          i=x)
#' })
#' @tests
#' data("ref_health_pancrea_cor_list")
#' querydata_list <- list(query_df)
#' test_query_predigs_score_list <- lapply(1:length(querydata_list),function(x){
#'     get_predigs_score_df(reference_cor_list = ref_health_pancrea_cor_list,
#'                                          querydata_cor_list = querydata_list,
#'                                          i=x)
#' })
#'expect_equal(dim(test_query_predigs_score_list[[1]]),dim(query_predigs_score_list[[1]]))
get_predigs_score_df <- function(reference_cor_list,querydata_cor_list,i){
  Tabula_cor_pre_df <- querydata_cor_list[[i]]
  reference_cor_pre_list <- reference_cor_list
  cor_matrix_list <- list()
  for (i in 1:length(reference_cor_pre_list)){
    reference_df <- reference_cor_pre_list[[i]]
    annotated_df <- Tabula_cor_pre_df
    intersect_gene <- intersect(colnames(reference_df),colnames(annotated_df))
    reference_df <- reference_df[,intersect_gene]
    annotated_df <- annotated_df[,intersect_gene]
    print(length(intersect_gene))
    cor_matrix <- cor(as.data.frame(t(reference_df)),as.data.frame(t(annotated_df)),method = "pearson")
    cor_matrix[is.na(cor_matrix)] <- 0
    cor_matrix_list[[i]] <- cor_matrix
  }
  cor_df <- c()

  for (i in 1:length(cor_matrix_list)){
    tmp <- reshape2::melt(cor_matrix_list[[i]])
    tmp$dataset <- paste0("dataset",i)
    cor_df <- rbind(cor_df,tmp)
  }

  colnames(cor_df) <- c("ReferenceCellType","QueryCellType","cor_value","RefDataset")
  cor_df$RefGroup <- cor_df$ReferenceCellType
  cor_df$RefGroup <- gsub("pancreas_","",cor_df$RefGroup)


  cor_integrated_list <- split.data.frame(cor_df,cor_df$QueryCellType)

  for (i in 1:length(cor_integrated_list)){
    tmp <- cor_integrated_list[[i]]
    tmp <- split.data.frame(tmp,tmp$RefGroup)
    tmp <- lapply(tmp,function(x){x <- as.vector(x$cor_value)})
    cor_integrated_list[[i]] <- tmp
  }

  cor_integrated_value <- cor_integrated_list
  ### calculate multi dataset cosine matrix
  for (i in 1:length(cor_integrated_list)){
    tmp <- cor_integrated_list[[i]]

    for (j in 1:length(tmp)){
      a <- cor_integrated_list[[i]][[j]]
      #print(a)
      b <- rep(1,length(a))
      c <- rep(0,length(a))
      c1 <- rbind(a,b)
      c2 <- rbind(b,c)
      value <-1-philentropy::distance(c1,method = "euclidean")/distance(c2,method = "euclidean")

      #print(value)
      cor_integrated_value[[i]][[j]] <- as.vector(value)
    }
  }

  cor_integrated_df <- as.data.frame(names(cor_integrated_value))

  for (i in 1:length(cor_integrated_value)){cor_integrated_df[i,2:(length(cor_integrated_value[[1]])+1)] <- unlist(cor_integrated_value[[i]])}
  colnames(cor_integrated_df) <- c("Celltype",names(cor_integrated_value[[1]]))
  rownames(cor_integrated_df) <- cor_integrated_df$Celltype

  #############   the end
  return(cor_integrated_df)
}


#' @title get_predigs_label_table()
#' @description  get predicted label table after predigs cell prediction
#' @details predicted cellular label output at cluster or single cell level
#' @param query_predigs_score_list a list of query-reference cell similarity score matrices after predigs prediction
#' @param clusterlevel TRUE or FALSE, which TRUE means cell type annotation at cluster-level ,FALSE means annotation at single cell-level
#'
#' @return
#' @export
#'
#' @examples
#' query_predigs_label_table_list <- get_predigs_label_table(query_predigs_score_list,clusterlevel = FALSE)
#' @tests
#' test <- get_predigs_label_table(query_predigs_score_list,clusterlevel = TRUE)
#' expect_equal(test[[1]],query_predigs_label_table_list[[1]])
get_predigs_label_table <- function(query_predigs_score_list,clusterlevel = TRUE){
  if (clusterlevel == TRUE) {
    query_predigs_label_table_list <- get_predigs_label_table_cluster(query_predigs_score_list)
  }

  if (clusterlevel == FALSE) {
    query_predigs_label_table_list <- get_predigs_label_table_sc(query_predigs_score_list)
  }

  ## the end
  return(query_predigs_label_table_list)
}



#' @title get_predigs_label_table_cluster
#' @description get label table of predigs grouped cell prediction
#' @description predicted cellular label output at cluster-level
#' @param query_predigs_score_list a list of query-reference cell similarity score matrices after predigs prediction
#'
#' @return
#' @export
#'
get_predigs_label_table_cluster <- function(query_predigs_score_list){
  query_predigs_score_list <- lapply(query_predigs_score_list,function(x){
    apply(x,1,function(t){colnames(x)[which.max(t)]})
  })
  for (i in 1:length(query_predigs_score_list)){
    new_df <- data.frame(
      group_name = as.character(names(query_predigs_score_list[[i]])),
      predigs_label = as.character(query_predigs_score_list[[i]])
    )
    query_predigs_score_list[[i]] <- new_df
  }
  return(query_predigs_score_list)
}

#' @title get_predigs_label_table_sc
#' @description get label table of predigs single cell prediction
#' @description predicted cellular label output at single cell-level
#' @param query_predigs_score_list a list of query-reference cell similarity score matrices after predigs prediction
#'
#' @return
#' @export
#'
get_predigs_label_table_sc <- function(query_predigs_score_list){
  query_predigs_score_list <- lapply(query_predigs_score_list,function(x){
    apply(x,1,function(t){colnames(x)[which.max(t)]})
  })
  return(query_predigs_score_list)
}

#' @title get_predigs_label_visualization
#' @description The visualization of query-reference cell similarity score matrix with predigs prediction
#' @param query_predigs_score_df the data frame of query-reference cell similarity score matrix with predigs prediction
#'
#' @export
#' @importFrom ComplexHeatmap Heatmap rowAnnotation HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.text
#' @examples
#' get_predigs_label_visualization(query_predigs_score_list[[1]])
get_predigs_label_visualization <- function(query_predigs_score_df){
  col_fun2 = circlize::colorRamp2(c(-1, 0.3, 0.5,0.7,0.9), c("blue", "white","orange","red","#660000"))
  cell_fun_all = function(j, i, x, y, width, height, fill){
    grid::grid.text(sprintf("%.3f", query_predigs_score_df[,2:length(colnames(query_predigs_score_df))][i,j]), x, y, gp = gpar(fontsize = 10))
  }
  row_ha = rowAnnotation(ReferenceCellType = rep("query",length(rownames(query_predigs_score_df))),
                         col = list(ReferenceCellType=c("query"="white")),
                         annotation_name_gp = gpar(fontface = "bold"),
                         show_legend = c(bar =FALSE))
  column_ha = HeatmapAnnotation(QueryCellGroup = rep("reference",length(colnames(query_predigs_score_df))-1),
                          col = list(QueryCellGroup=c("reference"="white")),
                          annotation_name_gp = gpar(fontface = "bold"),
                          show_legend = c(bar =FALSE))
  ComplexHeatmap::Heatmap(as.matrix(query_predigs_score_df[,2:length(colnames(query_predigs_score_df))]),
                          cell_fun = cell_fun_all,
                          col = col_fun2,
                          heatmap_legend_param =list(title = "Simiarlity Score"),
                          top_annotation = column_ha,
                          left_annotation = row_ha
                          )

}

#' @title get_predigs_label_visualization_sankey
#' @description The visualization of querygroup-reference cell label after predigs prediction
#' @param df the data frame of querygroup-reference cell label data frame after predigs prediction
#' @param GroupNumberBalance choose whether viualize the querygroup-reference releationship in a balanced query-group number ,with FALSE at default.
#' @param cols a vector of color palettes ,with NULL at default
#' @export
#' @import ggalluvial
#' @importFrom circlize colorRamp2
get_predigs_label_visualization_sankey <- function(df,GroupNumberBalance = FALSE,cols= NULL){
  # load data frame
  data <- df
  # parameter setting
  PredictedCellTypeCellType <- colnames(data)[1]
  CellGroupName <- colnames(data)[2]
  colnames(data) <- c("PredictedCellType","group")
  # preprocessing
  data$value <- 1
  data$blue <- paste0(data[,1],"|",data[,2])
  data <- aggregate.data.frame(data$value,by=list(data$blue),FUN=sum)
  data$group<- as.vector(unlist(lapply(data$Group.1,function(x){x <- strsplit(x,"[|]")[[1]][2]})))
  data$PredictedCellType <- as.vector(unlist(lapply(data$Group.1,function(x){x <- strsplit(x,"[|]")[[1]][1]})))
  colnames(data)[1:2] <- c("glue","CellNumber")
  data_sum <- aggregate(data$CellNumber,by=list(data[,"group"]),FUN=sum)
  data$CellProportion <- 0
  for (i in 1:length(rownames(data))){
    data$CellProportion[i] <-  data$CellNumber[i]/data_sum[which(data_sum$Group.1==data[i,"group"]),"x"]
  }
  # visualization setting
  if (!(is.null(cols))) {
    cols = cols
  }

  if (is.null(cols)) {
    cols <- sample(my36colors,length(levels(as.factor(data$group))))
  }


  # start visualization
  if (GroupNumberBalance == FALSE){
    ggplot(data = data,
           aes(axis2 = PredictedCellType, axis1 = group,
               y = CellNumber)) +
      scale_x_discrete(limits = c(CellGroupName,PredictedCellTypeCellType), expand = c(.2, .05)) +
      #xlab("Demographic") +
      geom_alluvium(aes(fill = group)) +
      geom_stratum(alpha = .0, width=1/3,color = "black") +
      scale_fill_brewer(type = "qual", palette = "Set1")+
      geom_text(stat = "stratum",
                aes(label = after_stat(stratum)))+
      #geom_flow(stat = "alluvium", lode.guidance = "zigzag") +
      theme_minimal() +
      ggtitle(paste0("The Sankey Diagram of Cell Assignment, by ",CellGroupName," and ",PredictedCellTypeCellType))+
      scale_fill_manual(values = cols)
  }

  if (GroupNumberBalance == TRUE){
    print("done!")
    ggplot(data = data,
           aes(axis1 = group, axis2 = PredictedCellType,
               y = CellProportion)) +
      scale_x_discrete(limits = c(CellGroupName,PredictedCellTypeCellType), expand = c(.2, .05)) +
      #xlab("Demographic") +
      geom_alluvium(aes(fill = group)) +
      geom_stratum(alpha = .0, width=1/3,color = "black") +
      scale_fill_brewer(type = "qual", palette = "Set1")+
      geom_text(stat = "stratum",
                aes(label = after_stat(stratum)))+
      #geom_flow(stat = "alluvium", lode.guidance = "zigzag") +
      theme_minimal() +
      ggtitle(paste0("The Sankey Diagram of Cell Assignment, by ",CellGroupName," and ",PredictedCellTypeCellType))+
      scale_fill_manual(values = cols)

  }
  # the end
}
