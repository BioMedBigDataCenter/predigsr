

#' @title predigs reference datasets list loading
#' @description load a list of multiple reference datasets
#' @details The reference list for sepcific organ with specific tissue types in the human digestive system
#' @param organtype specific organ type ,i.e.,intestine/pancrea/stomach/liver/esophagus
#' @param tissuetype specific tissue type ,i.e.,cancer/health
#'
#' @return
#' @export
#'
#' @examples
#' reference_dataset <- get_reference_setname("liver","health")
#' @tests
#' organtype = "liver"
#' tissuetype = "health"
#'
#' datasetname = get_reference_setname(organtype,tissuetype)
#' expect_equal(datasetname,"ref_health_liver_cor_list")
#'
get_reference_setname <- function(organtype=NULL,tissuetype=NULL){
  if(is.null(organtype)) {
    stop("Please provide one organ type (i.e.,intestine/pancrea/stomach/liver/esophagus)")
  }
  if(is.null(tissuetype)) {
    stop("Please provide one tissue type (i.e.,cancer/health)")
  }
  dataname <- paste0("ref_",tissuetype,"_",organtype,"_cor_list")
  print(dataname)
  return(dataname)
}
