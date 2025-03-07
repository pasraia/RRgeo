## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning=FALSE,
  message=FALSE
)

require(RRgeo)
require(rnaturalearth)
require(ggplot2)
require(terra)
require(sf)
options(rmarkdown.html_vignette.check_title = FALSE)
load(file.path(dirname(getwd()),"inst","exdata","Ursus_occurrences.RDa"))
rast(file.path(dirname(getwd()),"inst","exdata","U.arctos_suitability.tif"))->map1
rast(file.path(dirname(getwd()),"inst","exdata","U.maritimus_suitability.tif"))->map2

## ----loadoccs1, eval=FALSE----------------------------------------------------
#  library(RRgeo)
#  library(terra)
#  library(sf)
#  
#  rast(system.file("exdata/U.arctos_suitability.tif", package="RRgeo"))->map1
#  rast(system.file("exdata/U.maritimus_suitability.tif", package="RRgeo"))->map2
#  load(system.file("exdata/Ursus_occurrences.RDa", package="RRgeo"))

## ----loadoccs2----------------------------------------------------------------
list(Ursus_arctos=map1,Ursus_maritimus=map2)->pred
list(Ursus_arctos=occs_arctos,Ursus_maritimus=occs_marit)->occs
head(occs$Ursus_arctos)


## ----RRphylogeo,eval=FALSE----------------------------------------------------
#  setwd("YOUR_DIRECTORY")
#  
#  RRphylogeography(spec1="Ursus_arctos",
#                   spec2="Ursus_maritimus",
#                   pred=pred,
#                   occs=occs,
#                   aggr=2,
#                   time_col="TIME_factor",
#                   kde_inversion=FALSE,
#                   resistance_map=NULL,
#                   clust=0.5,
#                   plot=FALSE,
#                   mask_for_pred=NULL,
#                   standardize=TRUE,
#                   output.dir='.')
#  

## ----RPO,out.width='95%',dpi=300,echo=FALSE-----------------------------------
knitr::include_graphics("RPO_map.png")

