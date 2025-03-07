## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}

if (!requireNamespace("kableExtra", quietly = TRUE)) {
   warning(call. = FALSE, "kableExtra not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align='center',
  fig.width=6,
  fig.height=4,
  warning=FALSE,
  message=FALSE,
  out.width='95%',
  dpi=300
)

require(RRgeo)
require(rnaturalearth)
require(ggplot2)
require(ggtext)
require(terra)
require(sf)
require(ape)
require(viridis)

options(rmarkdown.html_vignette.check_title = FALSE)
load("enphylo_vignette.Rda")
# read.tree("Eucopdata_tree.txt")->tree
read.tree(file.path(dirname(getwd()),"inst","exdata","Eucopdata_tree.txt"))->tree
rast(file.path(dirname(getwd()),"inst","exdata","X35kya.tif"))->map35
rast(file.path(dirname(getwd()),"inst","exdata","Suit_Vulpes.tif"))->map1
rast(file.path(dirname(getwd()),"inst","exdata","Suit_Ursus.tif"))->map2
as.data.frame(map1,xy=T)->map1
as.data.frame(map2,xy=T)->map2

## ----loadstuff, eval=FALSE----------------------------------------------------
#  library(RRgeo)
#  library(terra)
#  library(sf)
#  library(ape)
#  
#  # datG<-lapply(grep(".gpkg",list.files(),value=TRUE),st_read)
#  # names(datG)<-sapply(strsplit(grep(".gpkg",list.files(),value=TRUE),"_"),"[[",1)
#  # dat<-lapply(datG,function(x) x[,c("OBS","age","bio1", "bio4", "bio11", "bio19")])
#  setwd("YOUR_DIRECTORY")
#  load(url("https://zenodo.org/records/14936297/files/dat.Rda?download=1"))
#  read.tree(system.file("exdata/Eucopdata_tree.txt", package="RRgeo"))->tree
#  tree$tip.label<-gsub("_"," ",tree$tip.label)

## ----dat----------------------------------------------------------------------
head(dat[1])

## ----enmod, eval=FALSE--------------------------------------------------------
#  rast(system.file("exdata/X35kya.tif", package="RRgeo"))->map35
#  project(map35,st_crs(dat[[1]])$proj4string,res = 50000)->map
#  
#  ENphylo_modeling(input_data=dat,
#                   tree=tree,
#                   input_mask=map[[1]],
#                   obs_col="OBS",
#                   time_col="age",
#                   min_occ_enfa=15,
#                   boot_test_perc=20,
#                   boot_reps=10,
#                   swap.args=list(nsim=10,si=0.2,si2=0.2),
#                   eval.args=list(eval_metric_for_imputation="AUC",
#                                  eval_threshold=0.7,
#                                  output_options="best"),
#                   clust=0.5,
#                   output.dir='.')
#  

## ----raewtgqa, eval=FALSE-----------------------------------------------------
#  library(rnaturalearth)
#  ne_countries(returnclass = "sf")->globalmap
#  subset(globalmap,continent=="North America")->ame_map
#  
#  map35[[c("bio1","bio4","bio11","bio19")]]->newmap
#  crop(newmap,ext(ame_map))->newmap
#  project(newmap,st_crs(dat[[1]])$proj4string,res = 50000)->newmap
#  
#  
#  getENphylo_results(input.dir ='.',
#                     mods="all",
#                     species_name=c("Vulpes velox","Ursus maritimus"))->mod
#  
#  ENphylo_prediction(object = mod,
#                     newdata = newmap,
#                     convert.to.suitability = TRUE,
#                     output.dir='.',
#                     proj_name="proj_example")
#  

## ----predmap, eval=FALSE------------------------------------------------------
#  library(ggplot2)
#  library(ggtext)
#  library(viridis)
#  
#  rast(system.file("exdata/Suit_Vulpes.tif", package="RRgeo"))->map1
#  rast(system.file("exdata/Suit_Ursus.tif", package="RRgeo"))->map2
#  as.data.frame(map1,xy=T)->map1
#  as.data.frame(map2,xy=T)->map2
#  
#  # rast("./ENphylo_prediction/Vulpes velox/proj_example/Suitability.tif")->map1
#  # rast("./ENphylo_prediction/Ursus maritimus/proj_example/Suitability.tif")->map2
#  

## ----predplot,fig.show='hold'-------------------------------------------------
p1<-ggplot(map1,aes(x=x,y=y,fill=Suitability_swap_9))+
  geom_tile()+
  scale_fill_viridis(name = "Suitability")+
  labs(title="*Vulpes velox* at 35 kya")+
  theme(panel.background = element_rect(fill="aliceblue",colour = "black"),
        panel.grid = element_blank(),
        axis.text= element_text(size=10),
        axis.title = element_blank(),
        plot.title = element_markdown(size=12,hjust=0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

p2<-ggplot(map2,aes(x=x,y=y,fill=Suitability))+
  geom_tile()+
  scale_fill_viridis()+
  labs(title="*Ursus maritimus* at 35 kya")+
  theme(panel.background = element_rect(fill="aliceblue",colour = "black"),
        panel.grid = element_blank(),
        axis.text=element_text(size=10),
        axis.title = element_blank(),
        plot.title = element_markdown(size=12,hjust=0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

plot(p1)
plot(p2)

