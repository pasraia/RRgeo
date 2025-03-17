#'@title Find species area of origin
#'@description The function integrates phylogenetic and geographical data (i.e.,
#'  habitat suitability maps), along with tools to model species' movements
#'  across landscapes. This integration allows for inference of the most
#'  probable area of origin (i.e. speciation) or regions of historical contact
#'  across time and space between a pair of target species.
#'@usage RRphylogeography(spec1,spec2,pred,occs,aggr=NULL,time_col=NULL,
#'  kde_inversion=FALSE,resistance_map=NULL,th=0.5,clust=0.5,plot=FALSE,
#'  mask_for_pred=NULL,standardize=TRUE,output.dir='.')
#'@param spec1,spec2 character. The names of the sister species whose area of
#'  origin should be inferred.
#'@param pred a list of two \code{SpatRaster} objects containing the prediction
#'  maps in logistic output generated through any Species Distribution Model
#'  technique. List names must correspond to \code{spec1} and \code{spec2}.
#'@param occs a list of two \code{sf::data.frame} objects containing species
#'  occurrence data only in binary format (exclusively ones for presence). List
#'  names must correspond to \code{spec1} and \code{spec2}.
#'@param aggr positive integer. Aggregation factor expressed as number of cells
#'  in each direction to be aggregated by averaging cell values (optional).
#'@param time_col character. Name of the \code{occs} column containing the time
#'  intervals associated to each species occurrence (optional).
#'@param kde_inversion logical. If \code{TRUE} and \code{time_col} is provided,
#'  kernel density is estimated by inverting the weights associated to the
#'  occurrences of the oldest species.
#'@param resistance_map a \code{SpatRaster} object to be optionally provided as
#'  a cost surface map to account for the difficulty of moving in each direction
#'  across the study area.
#'@param th numeric. The threshold value to define most suitable cells of the
#'  species pair as predicted via SDMs.
#'@param clust numeric. The proportion of the proportion of clusters to be used
#'  in parallel computing. Default is 0.5. If \code{NULL}, parallel computing is
#'  disabled.
#'@param plot logical. If \code{TRUE}, the area of the origin (RPO) maps with
#'  and without the kernel estimation factors are plotted.
#'@param mask_for_pred a \code{SpatRaster} object representing the geographical
#'  extent of the study area. This is used to crop the prediction maps define
#'  the area of interest for the analysis (optional).
#'@param standardize logical. If \code{TRUE}, the Relative Probability (RPO)
#'  values are standardized between 0 and 1.
#'@param output.dir character. The file path wherein \code{RRphylogeography}
#'  creates a new folder to store the outputs. This new folder is renamed by
#'  concatenating the names of the species pair.
#'@importFrom terra resample values merge mean
#'@importFrom sf st_as_sf st_length
#'@importFrom leastcostpath create_cs force_isotropy create_lcp
#'@importFrom doSNOW registerDoSNOW
#'@importFrom utils txtProgressBar setTxtProgressBar
#'@author Alessandro Mondanaro, Mirko Di Febbraro, Silvia Castiglione, Carmela
#'  Serio, Marina Melchionna, Pasquale Raia
#'@details \code{RRphylogeography} identifies the most suitable cells for both
#'  target species by relying on the \code{th} value. This threshold represents
#'  a numeric quantile so that any habitat suitability value greater than 0 that
#'  exceeds the value is considered to belong to the most suitable cells for the
#'  species. Conversely, the cells having habitat suitability values lower than
#'  \code{th} are excluded from distance calculation.
#'@seealso  \href{../doc/RRphylogeo.html}{\code{RRphylogeography} vignette}
#'@return A list of \code{SpatRaster} objects which includes the area of the
#'  origin (RPO) and both relative probability (RPO) maps for the species pair
#'  calculated for each layer in the prediction maps.
#'@export
#'@references Mondanaro, A., Castiglione, S., Di Febbraro, M., Timmermann, A.,
#'  Girardi, G., Melchionna, M., Serio, C., Belfiore, A.M., & Raia, P. (2025).
#'  RRphylogeography: A new method to find the area of origin of species and the
#'  history of past contacts between species. \emph{Methods in Ecology and
#'  Evolution}, 16: 546-557. 10.1111/2041-210X.14478
#'@examples
#' \dontrun{
#'
#' library(RRgeo)
#' library(terra)
#' library(sf)
#'
#' setwd("YOUR_DIRECTORY")
#' getwd()->main.dir
#'
#' rast(system.file("exdata/U.arctos_suitability.tif", package="RRgeo"))->map1
#' rast(system.file("exdata/U.maritimus_suitability.tif", package="RRgeo"))->map2
#' load(system.file("exdata/Ursus_occurrences.RDa"), package="RRgeo"))
#' list(Ursus_arctos=map1,Ursus_maritimus=map2)->pred
#' list(Ursus_arctos=occs_arctos,Ursus_maritimus=occs_marit)->occs
#'
#' RRphylogeography(spec1="Ursus_arctos",
#'                  spec2="Ursus_maritimus",
#'                  pred=pred,
#'                  occs=occs,
#'                  aggr=2,
#'                  time_col="TIME_factor",
#'                  kde_inversion=FALSE,
#'                  resistance_map=NULL,
#'                  clust=0.5,
#'                  plot=FALSE,
#'                  mask_for_pred=NULL,
#'                  standardize=TRUE,
#'                  output.dir=main.dir)
#'
#'}



RRphylogeography<-function(spec1,
                     spec2,
                     pred,
                     occs,
                     aggr=NULL,
                     time_col=NULL,
                     kde_inversion=FALSE,
                     resistance_map=NULL,
                     th=0.5,
                     clust=0.5,
                     plot=FALSE,
                     mask_for_pred=NULL,
                     standardize=TRUE,
                     output.dir='.'){

  # require(scales)
  # require(sf)
  # require(parallel)
  # require(spatialEco)
  # require(doParallel)
  # require(ggplot2)
  # require(cowplot)
  # require(doSNOW)
  # require(terra)
  # require(leastcostpath)

  if (plot&&(!requireNamespace("ggplot2", quietly = TRUE)|!requireNamespace("cowplot", quietly = TRUE))) {
    stop("Packages \"ggplot2\" and \"cowplot\" needed for plot=TRUE. Please install them.",
         call. = FALSE)
  }

  RRphylogeo_inner<-function(spec1,
                             spec2,
                             pred,
                             aggr=NULL,
                             occs,
                             time_col=NULL,
                             kde_inversion=FALSE,
                             resistance_map=NULL,
                             th=NULL,
                             clust=NULL,
                             plot=FALSE,
                             mask_for_pred=NULL,
                             standardize=TRUE,
                             out.dir=output.dir){

    dir.create(file.path(out.dir,paste(spec1,spec2,sep="_")),recursive=TRUE)

    pred[[which(names(pred)%in%spec1)]]->pred1
    pred[[which(names(pred)%in%spec2)]]->pred2

    occs[[which(names(occs)%in%spec1)]]->occs1
    occs[[which(names(occs)%in%spec2)]]->occs2

    terra::merge(pred1,pred2)->mm_mask
    if(!is.null(aggr)) terra::aggregate(mm_mask,aggr,mean)->mm_mask

    resample(pred1,mm_mask,method="bilinear")->suit1
    resample(pred2,mm_mask,method="bilinear")->suit2
    suit1[is.na(mm_mask)]<-NA
    suit2[is.na(mm_mask)]<-NA

    if (!is.null(mask_for_pred)){
      crop(suit1,ext(mask_for_pred))->suit1
      crop(suit2,ext(mask_for_pred))->suit2
    }

    #### KERNEL-DENSITY MAP####
    occs1[!is.na(extract(suit1,vect(occs1),ID=FALSE)),]->co1xx
    occs2[!is.na(extract(suit2,vect(occs2),ID=FALSE)),]->co2xx


    if(is.null(time_col)){

      suppressMessages({sf.kde.mod(x=co1xx,
                                   ref = suit1,
                                   standardize = TRUE)->k1})

      suppressMessages({sf.kde.mod(x=co2xx,
                                   ref = suit2,
                                   standardize = TRUE)->k2})

    } else {
      if(kde_inversion==TRUE){
        w<-which.max(c(max(as.numeric(as.data.frame(co1xx)[,grepl(time_col,names(co1xx))])),
                       max(as.numeric(as.data.frame(co2xx)[,grepl(time_col,names(co2xx))]))))

        y1<-scales::rescale(as.numeric(as.data.frame(co1xx)[,grepl(time_col,names(co1xx))]))
        y2<-scales::rescale(as.numeric(as.data.frame(co2xx)[,grepl(time_col,names(co2xx))]))
        if(w==1)y1<-(1-y1) else y2<-(1-y2)

        suppressMessages({sf.kde.mod(x=co1xx,
                                     y=y1,
                                     ref = suit1,
                                     standardize = TRUE)->k1})
        suppressMessages({sf.kde.mod(x=co2xx,
                                     y=y2,
                                     ref = suit2,
                                     standardize = TRUE)->k2})
      } else {

        suppressMessages({sf.kde.mod(x=co1xx,
                                     y=scales::rescale(as.numeric(as.data.frame(co1xx)[,grepl(time_col,names(co1xx))])),
                                     ref = suit1,
                                     standardize = TRUE)->k1})
        suppressMessages({sf.kde.mod(x=co2xx,
                                     y=scales::rescale(as.numeric(as.data.frame(co2xx)[,grepl(time_col,names(co2xx))])),
                                     ref = suit2,
                                     standardize = TRUE)->k2})

      }

    }

    mean(k1,k2)->kk


    #### SUITABILITY MAP ####
    suit1->sp_tot1
    suit2->sp_tot2

    if (max(values(sp_tot1),na.rm=TRUE)>0.1 & max(values(sp_tot2),na.rm=TRUE)>0.1){

      as.data.frame(sp_tot1,na.rm=TRUE)->nn1 ### as.data.frame from terra
      quantile(nn1[nn1>0.1],probs=th)->q1
      as.data.frame(sp_tot1,xy=TRUE,na.rm=TRUE)->sp_tot1 ### as.data.frame from terra
      sp_tot1[sp_tot1[[3]]>=q1,]->sp_tot1
      extract(suit1,sp_tot1[,c("x","y")],cells=TRUE)$cell->sp_tot1$cellID
      st_as_sf(sp_tot1,coords=c("x","y"))->sp_tot1
      suit1*(suit1[[1]]>=q1)->suit1_red

      as.data.frame(sp_tot2,na.rm=TRUE)->nn2 ### as.data.frame from terra
      quantile(nn2[nn2>0.1],probs=th)->q2
      as.data.frame(sp_tot2,xy=TRUE,na.rm=TRUE)->sp_tot2 ### as.data.frame from terra
      sp_tot2[sp_tot2[[3]]>=q2,]->sp_tot2
      extract(suit2,sp_tot2[,c("x","y")],cells=TRUE)$cell->sp_tot2$cellID
      st_as_sf(sp_tot2,coords=c("x","y"))->sp_tot2
      suit2*(suit2[[1]]>=q2)->suit2_red

      mean(suit1_red,suit2_red)->suit

      if (!is.null(mask_for_pred)) {crop(mm_mask,ext(mask_for_pred))->mm_mask}

      if (is.null(resistance_map)){
        create_cs(mm_mask,4)->tr_cor
        force_isotropy(tr_cor)->tr_cor

      } else {
        resistance_map<-resample(resistance_map, mm_mask)
        create_cs(mm_mask,4,dem=resistance_map)->tr_cor
        force_isotropy(tr_cor)->tr_cor
      }


      makeCluster(detectCores()*clust)->cl
      clusterEvalQ(cl,{
        library(leastcostpath)
        library(sf)
        library(terra)
      })

      prova3 <- pblapply(1:length(sp_tot1$cellID),
                         function(x){
                           print(x)
                           pp1 <- create_lcp(tr_cor, sp_tot1[x, ], sp_tot2,
                                             cost_distance = TRUE)
                           gd <- st_length(pp1)
                           kk <- data.frame(cellID1 = pp1$fromCell, cellID2 = pp1$toCell,
                                            distance = as.numeric(1/gd), cost = pp1$cost)

                           if(any(is.finite(kk$distance))) {

                             if (sp_tot1[x,]$cellID %in% sp_tot2$cellID) {
                               which(sp_tot1[x,]$cellID==sp_tot2$cellID)->q
                               kk[q,]$distance<-min(kk$distance[is.finite(kk$distance)],
                                                    na.rm=TRUE)
                               kk[q,]$cost<-min(kk$cost,na.rm=TRUE)
                             }

                             if(sum(kk$distance=="Inf")>0){
                               which(kk$distance=="Inf")->l
                               kk[l,]$distance<-max(kk$distance[is.finite(kk$distance)],
                                                    na.rm=TRUE)
                               kk[l,]$cost<-max(kk$cost,na.rm=TRUE)
                             }

                           } else kk<-NULL
                           kk
                         },cl=cl)
      stopCluster(cl)
      closeAllConnections()
      gc()

      do.call(rbind,prova3)->D_mat
      D_mat$distance[D_mat$distance==1e-51]<-0

      as.data.frame(suit1,xy=TRUE,cells=T)->grid1
      colnames(grid1)<-c("cell","x","y","layer")
      grid1$layer<-0
      stats::aggregate(distance~cellID1,FUN=mean,data=D_mat)->agg_sp1
      grid1$layer[grid1$cell%in%agg_sp1$cellID1]<-agg_sp1[,2]
      prova_d1<-suit1
      prova_d1[!is.na(prova_d1)]<-grid1$layer

      as.data.frame(suit2,xy=TRUE,cells=T)->grid2
      colnames(grid2)<-c("cell","x","y","layer")
      grid2$layer<-0
      stats::aggregate(distance~cellID2,FUN=mean,data=D_mat)->agg_sp2
      grid2$layer[grid2$cell%in%agg_sp2$cellID2]<-agg_sp2[,2]
      prova_d2<-suit2
      prova_d2[!is.na(prova_d2)]<-grid2$layer

      scales::rescale(values(prova_d1))->prova_d1[]
      scales::rescale(values(prova_d2))->prova_d2[]
      mean(prova_d1,prova_d2)->prova_dist

      RPO_combined=prova_dist*kk*suit
      RPO_combined_nok=prova_dist*suit

      RPO_sp1=suit1_red*k1*prova_d1
      RPO_sp1_nok=suit1_red*prova_d1

      RPO_sp2=suit2_red*k2*prova_d2
      RPO_sp2_nok=suit2_red*prova_d2

      TOT<-list(RPO_combined,
                RPO_combined_nok,
                RPO_sp1,
                RPO_sp1_nok,
                RPO_sp2,
                RPO_sp2_nok)
      names(TOT)<-c("RPO_combined",
                    "RPO_combined_nok",
                    "RPO_sp1",
                    "RPO_sp1_nok",
                    "RPO_sp2",
                    "RPO_sp2_nok")

      if (standardize){
        TOT<-lapply(TOT, function(x){

          x1<-scales::rescale(values(x))
          x[]<-x1
          x
        })
      }

      sp1<-c(suit1_red,k1,prova_d1,TOT$RPO_sp1,TOT$RPO_sp1_nok)
      sp2<-c(suit2_red,k2,prova_d2,TOT$RPO_sp2,TOT$RPO_sp2_nok)
      names(sp1)<-c("Suitability",
                    "Kernel_density",
                    "Proximity",
                    paste("RPO",spec1,sep="_"),
                    paste("RPO",spec1,"nok",sep="_"))
      names(sp2)<-c("Suitability",
                    "Kernel_density",
                    "Proximity",
                    paste("RPO",spec2,sep="_"),
                    paste("RPO",spec2,"nok",sep="_"))
      combined<-c(TOT$RPO_combined,TOT$RPO_combined_nok)
      names(combined)<-c("RPO_combined",
                         "RPO_combined_nok")

      list(name=names(suit1),
           RPO_sp1=sp1,
           RPO_sp2=sp2,
           RPO_combined=combined)->res_tot

      names(res_tot)[2:3]<-c(spec1,spec2)

      if (plot){

        pal<-c("white","#7FCDBB","#41B6C4","#2C7FB8","#00012E")

        as.data.frame(res_tot$RPO_combined,xy=TRUE)->dp1
        p1<-ggplot2::ggplot(co1xx)+
          ggplot2::geom_sf(colour="transparent",pch=21,fill="transparent",size=5)+
          ggplot2::geom_tile(data=dp1,ggplot2::aes(x=dp1$x,y=dp1$y,fill=dp1$RPO_combined))+
          ggplot2::scale_fill_gradientn(colours = pal)+
          ggplot2::theme(panel.background = ggplot2::element_rect(fill = "gray90",colour = "black"),
                         axis.title=ggplot2::element_blank(),
                         axis.text=ggplot2::element_text(size=13),
                         legend.position = "none")

        as.data.frame(res_tot$RPO_combined_nok,xy=TRUE)->dp2
        p2<-ggplot2::ggplot(co1xx)+
          ggplot2::geom_sf(colour="transparent",pch=21,fill="transparent",size=5)+
          ggplot2::geom_tile(data=dp2,ggplot2::aes(x=dp2$x,y=dp2$y,fill=dp2$RPO_combined_nok))+
          ggplot2::scale_fill_gradientn(colours = pal)+
          ggplot2::theme(panel.background = ggplot2::element_rect(fill = "gray90",colour = "black"),
                         axis.title=ggplot2::element_blank(),
                         axis.text=ggplot2::element_text(size=13),
                         legend.position = "none")

        cowplot::plot_grid(p1,p2,ncol=2,align="hv")->pp
        print(pp)
      }

      writeRaster(res_tot[[2]],
                  paste(out.dir,"/",paste(spec1,spec2,sep="_"),"/",
                        res_tot$name,"_",spec1,".tif",sep=""),overwrite=TRUE)

      writeRaster(res_tot[[3]],
                  paste(out.dir,"/",paste(spec1,spec2,sep="_"),"/",
                        res_tot$name,"_",spec2,".tif",sep=""),overwrite=TRUE)

      writeRaster(res_tot[[4]],
                  paste(out.dir,"/",paste(spec1,spec2,sep="_"),"/",
                        res_tot$name,"_combined.tif",sep=""),overwrite=TRUE)


    } else res_tot<-NULL


    return(res_tot)
  }


  outcome<-lapply(1:nlyr(pred[[1]]), function(layer){
    print(layer)
    RRphylogeo_inner(spec1,
                     spec2,
                     pred=lapply(pred, "[[", layer),
                     aggr=aggr,
                     occs=occs,
                     time_col=time_col,
                     th=th,
                     kde_inversion=kde_inversion,
                     resistance_map=resistance_map,
                     mask_for_pred=mask_for_pred,
                     clust=clust,
                     standardize=standardize,
                     plot=plot)
  })
  return(outcome)
}
