## ----include = FALSE----------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
     !rmarkdown::pandoc_available()) {
   warning(call. = FALSE, "Pandoc not found, the vignettes is not built")
   knitr::knit_exit()
}

misspacks<-sapply(c("rnaturalearth","ggplot2","viridis","ggtext"),requireNamespace,quietly=TRUE)
if(any(!misspacks)){
  warning(call. = FALSE,paste(names(misspacks)[which(!misspacks)],collapse=", "), "not found, the vignettes is not built")
   knitr::knit_exit()
}

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning=FALSE,
  message=FALSE,
  out.width='95%',
  fig.align='center',
  fig.width=6,
  fig.height=4,
  dpi=300,
  fig.pos='h'
)

require(RRgeo)
require(rnaturalearth)
require(ggplot2)
require(terra)
require(sf)
require(viridis)
require(ggtext)
options(rmarkdown.html_vignette.check_title = FALSE)
load(file.path(dirname(getwd()),"inst","exdata","Ursus_occs.Rda"))
load(file.path(dirname(getwd()),"inst","exdata","outputs.Rda"))
rast(file.path(dirname(getwd()),"inst","exdata","X35kya.tif"))->map35
map35[[1]]->map

## ----loadoccs1, eval=FALSE----------------------------------------------------
#  library(RRgeo)
#  library(terra)
#  library(sf)
#  library(rnaturalearth)
#  library(ggplot2)
#  library(ggtext)
#  library(viridis)
#  
#  
#  load(system.file("exdata/Ursus_occs.Rda", package="RRgeo"))

## ----loadoccs2----------------------------------------------------------------
ne_countries()->globalmap
subset(globalmap,continent=="Europe")->euromap
sf_use_s2(FALSE)
st_crop(euromap, xmin=-20, ymin=2 , xmax=50, ymax=75)->euromap
st_crop(occs, xmin=-20, ymin=2 , xmax=50, ymax=75)->occs
st_transform(euromap, st_crs("ESRI:54009"))->euromap
st_transform(occs, st_crs("ESRI:54009"))->occs

p1<-ggplot()+
  geom_sf(data=euromap,col="grey40",fill="white")+
  geom_sf(data=occs,fill="darkorchid2",size=2.5,color="black",pch=21)+
  labs(title="*Ursus spelaeus* occurrences")+
  theme(panel.background = element_rect(fill="aliceblue",colour = "black"),
        panel.grid = element_blank(),
        axis.text=element_text(size=10),
        plot.title = element_markdown(size=12,hjust=0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

plot(p1)


## ----mcp,out.width='47%',fig.width=4,fig.height=4,dpi=300,fig.show='hold',fig.align='default'----
buff<-0.1 # this is the element "buff" within the argument `bk_points`
sf::st_convex_hull(st_union(occs))->pol
max(st_distance(occs))*buff->buf 
st_buffer(pol,dist=as.numeric(buf))->pol_with_buffer

p2<-p1+
  geom_sf(data=pol,col="forestgreen",fill="transparent",lwd=1)+
  labs(title="*Ursus spelaeus* occurrences and MCP")+
  theme(panel.background = element_rect(fill="aliceblue",colour = "black"),
        panel.grid = element_blank(),
        axis.text= element_text(size=10),
        plot.title = element_markdown(size=12,hjust=0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
p3<-p1+
  geom_sf(data=pol,col="forestgreen",fill="transparent",lwd=1,linetype="dashed")+
  geom_sf(data=pol_with_buffer,col="forestgreen",fill="transparent",lwd=1)+
  labs(title="*Ursus spelaeus* occurrences and <br> MCP with buffer area")+
  theme(panel.background = element_rect(fill="aliceblue",colour = "black"),
        panel.grid = element_blank(),
        axis.text= element_text(size=10),
        plot.title = element_markdown(size=12,hjust=0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

plot(p2)
plot(p3)

## ----eval=FALSE---------------------------------------------------------------
#  curl::curl_download("https://zenodo.org/records/14998748/files/X35kya.tif?download=1",
#                     destfile = "X35kya.tif", quiet = FALSE)
#  rast("X35kya.tif")->map35
#  map35[[1]]->map

## ----density------------------------------------------------------------------
project(map,st_crs(pol_with_buffer)$proj4string,res = 50000)->map
mask(crop(map,vect(pol_with_buffer)),vect(pol_with_buffer))->map
RRgeo:::density_background(pres.locs=occs, MASK=map, rm.pres=TRUE)->dens.ras
as.data.frame(dens.ras,xy=TRUE)->dens.probs

p_dens<-ggplot()+
  geom_sf(data=euromap)+
  geom_tile(data=dens.probs,aes(x=x,y=y,fill=lyr.1))+
  scale_fill_viridis()+
  geom_sf(data=pol_with_buffer,col="forestgreen",fill="transparent",lwd=1)+
  labs(title="*Ursus spelaeus* density map",
       fill = "Probability")+
  theme(panel.background = element_rect(fill="aliceblue",colour = "black"),
        panel.grid = element_blank(),
        axis.text= element_text(size=10),
        axis.title = element_blank(),
        plot.title = element_markdown(size=12,hjust=0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
plot(p_dens)

## ----loadpseudoabsence, eval=FALSE--------------------------------------------
#  load(system.file("exdata/outputs.Rda",package="RRgeo"))

## ----pseudoabsence, echo=-1---------------------------------------------------
load(file.path(dirname(getwd()),"inst","exdata","outputs.Rda"))
ggplot()+
  geom_sf(data=euromap,col="grey40",fill="white")+
  geom_sf(data=psabs.points,aes(fill="pseudoabsence"),size=2.5,color="black",pch=21)+
  geom_sf(data=pres.points,aes(fill="presence"),size=2.5,color="black",pch=21)+
  scale_fill_manual(values=c("darkorchid2","gold2"))+
  geom_sf(data=pol_with_buffer,col="forestgreen",fill="transparent",lwd=1)+
  labs(title="*Ursus spelaeus* at 35 kya")+
  theme(panel.background = element_rect(fill="aliceblue",colour = "black"),
        panel.grid = element_blank(),
        axis.text=element_text(size=10),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white", colour = "white"),
        plot.title = element_markdown(size=12,hjust=0.5),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

## ----eucop, eval=FALSE--------------------------------------------------------
#  species_vec<-c( "Canis latrans","Canis lupus","Gulo gulo","Lutra lutra",
#                  "Martes americana","Meles meles","Mustela erminea",
#                  "Mustela nivalis","Procyon lotor","Ursus arctos","Ursus ingressus",
#                  "Ursus maritimus","Ursus spelaeus","Vulpes velox","Vulpes vulpes" )
#  
#  setwd("YOUR_DIRECTORY")
#  eucop_data_preparation(input.dir=getwd(), species_name=species_vec, variables="bio",
#                         calibration=TRUE, combine.ages="mean", bk_points=list(),
#                         output.dir=getwd())

