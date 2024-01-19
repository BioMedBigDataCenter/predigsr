
<!-- README.md is generated from README.Rmd. Please edit that file -->
PreDigsR
=======

<img src="https://github.com/jymeng-monica/predigs/blob/master/inst/logo.r.png" width="200" height="200" align="right"/>


<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/hexSticker?color=green)](https://cran.r-project.org/package=hexSticker)
<!-- badges: end -->

PreDigs is an automated method for cell type assignment of human digestive system with single cell RNA sequencing data

* Species:  support both `human` and `mouse`.

* Level: support query with either `single cells` or known `cell group labels`.

## Installation

This is an [R](https://www.r-project.org/) package.

You can install the released version of predigs from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("predigs")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jymeng-monica/predigs")
```

## Library

``` r
suppressMessages(library(Seurat))
suppressMessages(library(harmony)) # Version 0.1.0 is recommended
suppressMessages(library(dplyr))
suppressMessages(library(predigs))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggalluvial)) # optional
suppressMessages(library(circlize))
```

## Load Query and Reference Data

### Load Query Seurat Object

The required query data format of predigs recommends `normalized` Seurat object,i.e., with [seurat@assays$RNA@data](https://swbioinf.github.io/scRNAseqInR_Doco/seuratobject.html).
```r
load("./NormalizedSeuratObject.RData") # querydata
```

The seurat object with only `count` matrix are supported to run Seurat Workflow:
``` r
load("./RawSeuratObject.RData") # querydata
querydata <- get_seurat_workflow(querydata,stim = "batch_group_label")
```

###   Preprocess Query Data before Predigs Prediction

The input includes `querydata` and the colnames of `cell_group_label` in `querydata@meta.data`(optional)
```r
# 1.For human query data at single cell levels
query_df <- get_query_preprocess(querydata,species = "human",clusterlevel = FALSE)

# 2.For human query data at group levels
query_df <- get_query_preprocess(querydata,species = "human",cellgroup="cell_ontology_class",clusterlevel = TRUE)

# 3.For mouse query data at single cell levels
query_df <- get_query_preprocess(querydata,species = "mouse",clusterlevel = FALSE)

# 4.For mouse query data at group levels
query_df <- get_query_preprocess(querydata,species = "mouse",cellgroup="cell_ontology_class",clusterlevel = TRUE)

```
User should generate both a `query_df list` and a `querydata list` who have multiple query datasets
```r
querydata_list <_ list(querydata1,querydata2)
query_df_list <- list(query_df_1,query_df_2)
```

### Load Preprocessed Reference data 

User are required to specify the `Organ Type` name and `Tissue Type` name for their scRNAseq query data.

* Organ Type : intestine/pancrea/stomach/liver/esophagus

* Tissue Type : cancer/health

User can obtain corresponding `internal reference dataset name` of predigs and load it in R
```r
reference_dataset_name <- get_reference_setname("pancrea","health")
> "ref_health_pancrea_cor_list"
data(ref_health_pancrea_cor_list,package="predigs")
```
## Cell Type Assignment with PreDigs

### Cell Similarity Score Calculation 

```r
 query_predigs_score_list <- lapply(1:length(querydata_list),function(x){
     get_predigs_score_df(reference_cor_list = ref_health_pancrea_cor_list,
                                          querydata_cor_list = querydata_list,i=x)
```

### Ultimate Predicted Cell Type Label Output

```r
# For data query at single cell level
 query_predigs_label_table_list <- get_predigs_label_table(query_predigs_score_list,clusterlevel = FALSE)
 
# For data query at group level
 query_predigs_label_table_list <- get_predigs_label_table(query_predigs_score_list,clusterlevel = TRUE)
```

### The visualization of the Assignment Result

* Heatmap Plot
```r
get_predigs_label_visualization(query_predigs_score_list[[1]])
```
<img src="https://github.com/jymeng-monica/predigs/blob/master/inst/heatmapExm.png" width="100%" />


* Sankey Digram Plot
In this process,user are required to prepare a two-column data frame with the first column of `predicted cell label` and the second column of their `interested group` such as `ground truth` of cells. And the length of `rownames` should be the same as `cell numbers` of query data.

```r
head(sankry_df)
#>          predicted_type           ground_truth
#> 1      Acinar.cell pancreatic acinar cell
#> 2        CD8T.cell                 t cell
#> 3 Endothelial.cell       endothelial cell


# plot 
get_predigs_label_visualization_sankey(sankey_df,GroupNumberBalance = FALSE,cols= NULL)
```
<img src="https://github.com/jymeng-monica/predigs/blob/master/inst/sankey1.png" width="100%" />

User can browse the Sankey plot with a balanced group number like this:
```r
get_predigs_label_visualization_sankey(sankey_df,GroupNumberBalance = TRUE,cols= NULL)
```
<img src="https://github.com/jymeng-monica/predigs/blob/master/inst/sankey2.png" width="100%" />

User can specify `cols` if they need by :
```r
cols <- c("#E5D2DD","#53A85F","#F1BB72","#F3B1A0","#D6E7A3"....."#57C3F3")
#get_predigs_label_visualization_sankey(sankey_df,GroupNumberBalance = TRUE,cols= cols)
#get_predigs_label_visualization_sankey(sankey_df,GroupNumberBalance = FALSE,cols= NULL)
```

_____________________________________________________________
&copy; 2022.12.31 Monica Meng.Zhang Lab. All rights reserved.


