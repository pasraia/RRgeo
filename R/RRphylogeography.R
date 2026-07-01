#'@title Find species area of origin
#'@description The function integrates phylogenetic and geographical data (i.e.,
#'  habitat suitability maps), along with tools to model species' movements
#'  across landscapes. This integration allows for inference of the most
#'  probable area of origin (i.e. speciation) or regions of historical contact
#'  across time and space between a pair of target species.
#'@usage
#'RRphylogeography(spec1,spec2,pred,occs,aggr=NULL,time_col=NULL,weights=c(0.5,0.5),
#'kde_inversion=FALSE,resistance_map=NULL,th=0.5,clust=0.5,plot=FALSE,
#'mask_for_pred=NULL,standardize=TRUE,output.dir)
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
#'@param weights weights to account for the arithmetic (first value) and
#'  geometric (second value) means at calculating averaged suitability and
#'  estimate kernel density (see details).
#'@param kde_inversion logical. If \code{TRUE} and \code{time_col} is provided,
#'  kernel density is estimated by inverting the weights associated to the
#'  occurrences of the oldest species.
#'@param resistance_map an optional \code{SpatRaster} object representing a
#'  conductance matrix that numerically quantifies the resistance to move across
#'  a surface (0 indicating maximum resistance, 1 indicating minimum
#'  resistance).
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
#'@importFrom sf st_as_sf st_length st_linestring
#'@importFrom leastcostpath create_cs force_isotropy
#'@importFrom doSNOW registerDoSNOW
#'@importFrom utils txtProgressBar setTxtProgressBar
#'@author Alessandro Mondanaro, Mirko Di Febbraro, Silvia Castiglione, Carmela
#'  Serio, Marina Melchionna, Pasquale Raia
#'@details \code{RRphylogeography} identifies the most suitable cells for both
#'  target species by relying on the \code{th} value. This threshold represents
#'  a numeric quantile so that any habitat suitability value greater than 0 that
#'  exceeds the value is considered to belong to the most suitable cells for the
#'  species. Conversely, the cells having habitat suitability values lower than
#'  \code{th} are excluded from distance calculation. When averaging the habitat
#'  suitabilities and kernel densities of target species both the arithmetic
#'  mean and the geometric mean are computed. The final combined surfaces are
#'  defined as a weighted average of the two means, with weights summing to 1,
#'  according to the formula: weights [1]\*arithmetic mean +
#'  weights[2]\*geometric mean
#'@seealso  \href{../doc/RRphylogeo.html}{\code{RRphylogeography} vignette}
#'@return A list of \code{SpatRaster} objects which includes the area of the
#'  origin (AOO) and both relative probability (RPO) maps for the species pair
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
#' newwd<-tempdir()
#' # newwd<-"YOUR_DIRECTORY"
#'
#' rast(system.file("exdata/U.arctos_suitability.tif", package="RRgeo"))->map1
#' rast(system.file("exdata/U.maritimus_suitability.tif", package="RRgeo"))->map2
#' load(system.file("exdata/Ursus_occurrences.Rda", package="RRgeo"))
#' list(Ursus_arctos=map1,Ursus_maritimus=map2)->pred
#' list(Ursus_arctos=occs_arctos,Ursus_maritimus=occs_marit)->occs
#'
#' RRphylogeography(spec1="Ursus_arctos",
#'                  spec2="Ursus_maritimus",
#'                  pred=pred,
#'                  occs=occs,
#'                  aggr=5,
#'                  time_col="TIME_factor",
#'                  weights=c(0.5,0.5),
#'                  kde_inversion=FALSE,
#'                  resistance_map=NULL,
#'                  clust=NULL,
#'                  plot=FALSE,
#'                  mask_for_pred=NULL,
#'                  th=0.7,
#'                  standardize=TRUE,
#'                  output.dir=newwd)
#'}



RRphylogeography<-function(spec1,
                           spec2,
                           pred,
                           occs,
                           aggr=NULL,
                           time_col=NULL,
                           weights=c(0.5,0.5),
                           kde_inversion=FALSE,
                           resistance_map=NULL,
                           th=0.5,
                           clust=0.5,
                           plot=FALSE,
                           mask_for_pred=NULL,
                           standardize=TRUE,
                           output.dir){
  if (plot && (!requireNamespace("ggplot2", quietly = TRUE) |
               !requireNamespace("cowplot", quietly = TRUE)|
               !requireNamespace("igraph", quietly = TRUE))) {
    stop("Packages \"igraph\", \"ggplot2\" and \"cowplot\" needed for plot=TRUE. Please install them.",
         call. = FALSE)
  }

  if (is.null(output.dir))
    stop("argument \"output.dir\" is missing, with no default")
  if(!is.null(mask_for_pred)){
    if (!all(crs(pred[[1]],proj=TRUE)==crs(mask_for_pred,proj=TRUE),
             crs(pred[[2]],proj=TRUE)==crs(mask_for_pred,proj=TRUE))) {
      stop("CRS for the objects pred and mask_for_pred objects do not match.
         Please ensure all SpatRasters have the same projection and spatial resolution.")
    }
  }

  if (!is.null(resistance_map)){
    if (!all(crs(pred[[1]],proj=TRUE)==crs(resistance_map,proj=TRUE),
             crs(pred[[2]],proj=TRUE)==crs(resistance_map,proj=TRUE),
             all(res(pred[[1]])==res(resistance_map)),
             all(res(pred[[2]])==res(resistance_map)))) {
      stop("CRS and/or resolution for the objects pred and resistance_map do not match.
           Please ensure all SpatRasters have the same projection and spatial resolution.")
    }
  }

  RRphylogeo_inner <- function(spec1, spec2, pred, aggr = NULL,
                               occs, time_col = NULL, kde_inversion = FALSE, resistance_map = NULL,
                               th = NULL, clust = NULL, plot = FALSE, mask_for_pred = NULL,
                               standardize = TRUE, out.dir = output.dir) {
    dir.create(file.path(out.dir, paste(spec1, spec2, sep = "_")),
               recursive = TRUE)
    pred1 <- pred[[which(names(pred) %in% spec1)]]
    pred2 <- pred[[which(names(pred) %in% spec2)]]
    occs1 <- occs[[which(names(occs) %in% spec1)]]
    occs2 <- occs[[which(names(occs) %in% spec2)]]
    mm_mask <- terra::merge(pred1, pred2)
    if (!is.null(aggr))
      mm_mask <- terra::aggregate(mm_mask, aggr, mean,na.rm=TRUE)
    suit1 <- resample(pred1, mm_mask, method = "bilinear")
    suit2 <- resample(pred2, mm_mask, method = "bilinear")
    suit1[is.na(mm_mask)] <- NA
    suit2[is.na(mm_mask)] <- NA
    if (!is.null(mask_for_pred)) {
      suit1 <- crop(suit1, ext(mask_for_pred))
      suit2 <- crop(suit2, ext(mask_for_pred))
    }


    sp_tot1 <- suit1
    sp_tot2 <- suit2

    if (max(values(sp_tot1), na.rm = TRUE) > 0.1 & max(values(sp_tot2),
                                                       na.rm = TRUE) > 0.1) {
      nn1 <- as.data.frame(sp_tot1, na.rm = TRUE)
      q1 <- quantile(nn1[nn1 > 0.1], probs = th)
      sp_tot1 <- as.data.frame(sp_tot1, xy = TRUE, na.rm = TRUE)
      sp_tot1 <- sp_tot1[sp_tot1[[3]] >= q1, ]
      sp_tot1$cellID <- extract(suit1, sp_tot1[, c("x",
                                                   "y")], cells = TRUE)$cell
      sp_tot1 <- st_as_sf(sp_tot1, coords = c("x", "y"))
      suit1_red <- suit1 * (suit1[[1]] >= q1)
      nn2 <- as.data.frame(sp_tot2, na.rm = TRUE)
      q2 <- quantile(nn2[nn2 > 0.1], probs = th)
      sp_tot2 <- as.data.frame(sp_tot2, xy = TRUE, na.rm = TRUE)
      sp_tot2 <- sp_tot2[sp_tot2[[3]] >= q2, ]
      sp_tot2$cellID <- extract(suit2, sp_tot2[, c("x",
                                                   "y")], cells = TRUE)$cell
      sp_tot2 <- st_as_sf(sp_tot2, coords = c("x", "y"))
      suit2_red <- suit2 * (suit2[[1]] >= q2)
      mm_arith <- (suit1_red + suit2_red)/2
      mm_geom <- sqrt(suit1_red * suit2_red)
      suit <- weights[1] * mm_arith + weights[2] * mm_geom

      ## KDE ##
      co1xx <- occs1[!is.na(extract(suit1, vect(occs1), ID = FALSE)),
      ]
      co2xx <- occs2[!is.na(extract(suit2, vect(occs2), ID = FALSE)),
      ]
      if (is.null(time_col)) {
        suppressMessages({
          k1 <- sf.kde.mod(x = co1xx, ref = suit1, standardize = TRUE)
        })
        suppressMessages({
          k2 <- sf.kde.mod(x = co2xx, ref = suit2, standardize = TRUE)
        })
      }else {
        if (kde_inversion == TRUE) {
          w <- which.max(c(max(as.numeric(as.data.frame(co1xx)[,
                                                               grepl(time_col, names(co1xx))])), max(as.numeric(as.data.frame(co2xx)[,
                                                                                                                                     grepl(time_col, names(co2xx))]))))
          y1 <- scales::rescale(as.numeric(as.data.frame(co1xx)[,
                                                                grepl(time_col, names(co1xx))]))
          y2 <- scales::rescale(as.numeric(as.data.frame(co2xx)[,
                                                                grepl(time_col, names(co2xx))]))
          if (w == 1)
            y1 <- (1 - y1) else y2 <- (1 - y2)
          suppressMessages({
            k1 <- sf.kde.mod(x = co1xx, y = y1, ref = suit1,
                             standardize = TRUE)
          })
          suppressMessages({
            k2 <- sf.kde.mod(x = co2xx, y = y2, ref = suit2,
                             standardize = TRUE)
          })
        }else {
          suppressMessages({
            k1 <- sf.kde.mod(x = co1xx, y = scales::rescale(as.numeric(as.data.frame(co1xx)[,
                                                                                            grepl(time_col, names(co1xx))])), ref = suit1,
                             standardize = TRUE)
          })
          suppressMessages({
            k2 <- sf.kde.mod(x = co2xx, y = scales::rescale(as.numeric(as.data.frame(co2xx)[,
                                                                                            grepl(time_col, names(co2xx))])), ref = suit2,
                             standardize = TRUE)
          })
        }
      }

      k_arith <- (k1 + k2)/2
      k_geom <- sqrt(k1 * k2)
      kk <- weights[1] * k_arith + weights[2] * k_geom
      k1 <- mask(crop(k1, suit), suit)
      k2 <- mask(crop(k2, suit), suit)
      kk <- mask(crop(kk, suit), suit)
      if (!is.null(mask_for_pred)) {
        mm_mask <- crop(mm_mask, ext(mask_for_pred))

      }
      mm_mask[!is.na(mm_mask),]<-1
      if (is.null(resistance_map)) {
        tr_cor <- create_cs(mm_mask, 4)
        tr_cor <- force_isotropy(tr_cor)
      }else {
        resistance_map <- resample(resistance_map, mm_mask)
        mm_mask <- mm_mask*resistance_map
        tr_cor <- create_cs(mm_mask, 4)
        tr_cor <- force_isotropy(tr_cor)
      }
      cs_rast <- rast(nrow = tr_cor$nrow, ncol = tr_cor$ncol,
                      xmin = tr_cor$extent[1], xmax = tr_cor$extent[2],
                      ymin = tr_cor$extent[3], ymax = tr_cor$extent[4],
                      crs  = tr_cor$crs)
      cm_graph <- igraph::graph_from_adjacency_matrix(
        tr_cor$conductanceMatrix, mode = "directed", weighted = TRUE
      )
      igraph::E(cm_graph)$weight <- 1 / igraph::E(cm_graph)$weight

      if (!is.null(clust)) {
        cl <- makeCluster(detectCores() * clust)
        clusterEvalQ(cl, {
          library(leastcostpath)
          library(sf)
          library(terra)
        })
        prova3 <- pblapply(1:length(sp_tot1$cellID),
                           function(x) {
                             if(all(sp_tot2$cellID%in%sp_tot1[x, ])){
                               kk <- data.frame(cellID1=sp_tot1[x,]$cellID,
                                                cellID2=sp_tot2$cellID,
                                                distance=0,
                                                cost=1)
                             } else {
                               # NUOVA VERSIONE 23-01 ----------------------------------------------------
                               to_cell <- sp_tot2$cellID[sp_tot2$cellID != sp_tot1[x, ]$cellID]
                               from_cell <- sp_tot1[x, ]$cellID

                               sp <-  suppressMessages(igraph::shortest_paths(
                                 cm_graph, from = from_cell, to = to_cell,
                                 mode = "out", algorithm = "dijkstra"
                               ))

                               cost_vec <- as.numeric(igraph::distances(
                                 cm_graph, v = from_cell, to = to_cell,
                                 mode = "out", algorithm = "dijkstra"
                               ))

                               lcps_geom <- lapply(sp$vpath, function(v) {
                                 v <- as.integer(v)
                                 if (length(v)<= 1) {
                                   cell<-if(length(v)==0) from_cell else v[1]
                                   xy<-xyFromCell(cs_rast, cell)
                                   return(st_linestring(rbind(xy, xy)))
                                 }
                                 xy<-xyFromCell(cs_rast,v)
                                 st_linestring(xy)
                               })

                               lcps_sfc<-st_sfc(lcps_geom,crs=crs(cs_rast))
                               gd_vec<-as.numeric(st_length(lcps_sfc))

                               kk <- data.frame(cellID1 = from_cell,
                                                cellID2 = to_cell, distance = as.numeric(1/gd_vec),
                                                cost = cost_vec)
                               # FINE NUOVA VERSIONE 23-01 -----------------------------------------------
                               # sp_tot2[!sp_tot2$cellID%in%sp_tot1[x, ],]->sp_tot
                               # pp1 <- suppressMessages(create_lcp(tr_cor,
                               #                                    sp_tot1[x, ], sp_tot, cost_distance = TRUE))
                               # gd <- st_length(pp1)
                               # kk <- data.frame(cellID1 = pp1$fromCell,
                               #                  cellID2 = pp1$toCell, distance = as.numeric(1/gd),
                               #                  cost = pp1$cost)
                               if (any(is.finite(kk$distance))) {
                                 if (sp_tot1[x, ]$cellID %in% sp_tot2$cellID) {
                                   q <- which(sp_tot1[x, ]$cellID == sp_tot2$cellID)
                                   ss<-data.frame(cellID1=sp_tot1[x, ]$cellID,
                                                  cellID2=sp_tot2[q, ]$cellID,
                                                  distance=min(kk$distance[is.finite(kk$distance)],
                                                               na.rm = TRUE),
                                                  cost=min(kk$cost, na.rm = TRUE))
                                   rbind(kk,ss)->kk
                                 }
                                 if (sum(kk$distance == "Inf") > 0) {
                                   l <- which(kk$distance == "Inf")
                                   kk[l, ]$distance <- max(kk$distance[is.finite(kk$distance)],
                                                           na.rm = TRUE)
                                   kk[l, ]$cost <- max(kk$cost, na.rm = TRUE)
                                 }
                               }
                               else kk <- NULL
                             }
                             kk
                           }, cl = cl)
        stopCluster(cl)
        closeAllConnections()
        gc()
      }else {
        prova3 <- pblapply(1:length(sp_tot1$cellID),
                           function(x) {
                             if(all(sp_tot2$cellID%in%sp_tot1[x, ])){
                               kk <- data.frame(cellID1=sp_tot1[x,]$cellID,
                                                cellID2=sp_tot2$cellID,
                                                distance=0,
                                                cost=1)
                             } else {
                               # NUOVA VERSIONE 23-01 ----------------------------------------------------
                               to_cell <- sp_tot2$cellID[sp_tot2$cellID != sp_tot1[x, ]$cellID]
                               from_cell <- sp_tot1[x, ]$cellID

                               sp <-  suppressMessages(igraph::shortest_paths(
                                 cm_graph, from = from_cell, to = to_cell,
                                 mode = "out", algorithm = "dijkstra"
                               ))

                               cost_vec <- as.numeric(igraph::distances(
                                 cm_graph, v = from_cell, to = to_cell,
                                 mode = "out", algorithm = "dijkstra"
                               ))

                               lcps_geom <- lapply(sp$vpath, function(v) {
                                 v <- as.integer(v)
                                 if (length(v)<= 1) {
                                   cell<-if(length(v)==0) from_cell else v[1]
                                   xy<-xyFromCell(cs_rast, cell)
                                   return(st_linestring(rbind(xy, xy)))
                                 }
                                 xy<-xyFromCell(cs_rast,v)
                                 st_linestring(xy)
                               })

                               lcps_sfc<-st_sfc(lcps_geom,crs=crs(cs_rast))
                               gd_vec<-as.numeric(st_length(lcps_sfc))

                               kk <- data.frame(cellID1 = from_cell,
                                                cellID2 = to_cell, distance = as.numeric(1/gd_vec),
                                                cost = cost_vec)

                               # FINE NUOVA VERSIONE 23-01 ------------------------------------------
                               # sp_tot2[!sp_tot2$cellID%in%sp_tot1[x, ],]->sp_tot
                               # pp1 <- suppressMessages(create_lcp(tr_cor,
                               #                                    sp_tot1[x, ], sp_tot, cost_distance = TRUE))
                               # gd <- st_length(pp1)
                               # kk <- data.frame(cellID1 = pp1$fromCell,
                               #                  cellID2 = pp1$toCell, distance = as.numeric(1/gd),
                               #                  cost = pp1$cost)
                               if (any(is.finite(kk$distance))) {
                                 if (sp_tot1[x, ]$cellID %in% sp_tot2$cellID) {
                                   q <- which(sp_tot1[x, ]$cellID == sp_tot2$cellID)
                                   ss<-data.frame(cellID1=sp_tot1[x, ]$cellID,
                                                  cellID2=sp_tot2[q, ]$cellID,
                                                  distance=min(kk$distance[is.finite(kk$distance)],
                                                               na.rm = TRUE),
                                                  cost=min(kk$cost, na.rm = TRUE))
                                   rbind(kk,ss)->kk
                                 }
                                 if (sum(kk$distance == "Inf") > 0) {
                                   l <- which(kk$distance == "Inf")
                                   kk[l, ]$distance <- max(kk$distance[is.finite(kk$distance)],
                                                           na.rm = TRUE)
                                   kk[l, ]$cost <- max(kk$cost, na.rm = TRUE)
                                 }
                               }
                               else kk <- NULL
                             }
                             kk
                           })
      }

      if(!all(sapply(prova3,is.null))){
        D_mat <- do.call(rbind, prova3)
        D_mat$distance[D_mat$distance == 1e-51] <- 0
        grid1 <- as.data.frame(suit1, xy = TRUE, cells = T)
        colnames(grid1) <- c("cell", "x", "y", "layer")
        grid1$layer <- 0
        agg_sp1 <- stats::aggregate(distance ~ cellID1, FUN = mean,
                                    data = D_mat)
        grid1$layer[grid1$cell %in% agg_sp1$cellID1] <- agg_sp1[,
                                                                2]
        prova_d1 <- suit1
        prova_d1[!is.na(prova_d1)] <- grid1$layer
        grid2 <- as.data.frame(suit2, xy = TRUE, cells = T)
        colnames(grid2) <- c("cell", "x", "y", "layer")
        grid2$layer <- 0
        agg_sp2 <- stats::aggregate(distance ~ cellID2, FUN = mean,
                                    data = D_mat)
        grid2$layer[grid2$cell %in% agg_sp2$cellID2] <- agg_sp2[,
                                                                2]
        prova_d2 <- suit2
        prova_d2[!is.na(prova_d2)] <- grid2$layer
        prova_d1[] <- scales::rescale(values(prova_d1))
        prova_d2[] <- scales::rescale(values(prova_d2))
        prova_dist <- mean(prova_d1, prova_d2)
        RPO_combined = prova_dist * kk * suit
        RPO_combined_nok = prova_dist * suit
        RPO_sp1 = suit1_red * k1 * prova_d1
        RPO_sp1_nok = suit1_red * prova_d1
        RPO_sp2 = suit2_red * k2 * prova_d2
        RPO_sp2_nok = suit2_red * prova_d2
        TOT <- list(RPO_combined, RPO_combined_nok, RPO_sp1,
                    RPO_sp1_nok, RPO_sp2, RPO_sp2_nok)
        names(TOT) <- c("RPO_combined", "RPO_combined_nok",
                        "RPO_sp1", "RPO_sp1_nok", "RPO_sp2", "RPO_sp2_nok")
        if (standardize) {
          TOT <- lapply(TOT, function(x) {
            x1 <- scales::rescale(values(x))
            x[] <- x1
            x
          })
        }
        sp1 <- c(suit1_red, k1, prova_d1, TOT$RPO_sp1, TOT$RPO_sp1_nok)
        sp2 <- c(suit2_red, k2, prova_d2, TOT$RPO_sp2, TOT$RPO_sp2_nok)
        names(sp1) <- c("Suitability", "Kernel_density",
                        "Proximity", paste("RPO", spec1, sep = "_"),
                        paste("RPO", spec1, "nok", sep = "_"))
        names(sp2) <- c("Suitability", "Kernel_density",
                        "Proximity", paste("RPO", spec2, sep = "_"),
                        paste("RPO", spec2, "nok", sep = "_"))
        combined <- c(TOT$RPO_combined, TOT$RPO_combined_nok)
        names(combined) <- c("RPO_combined", "RPO_combined_nok")
        res_tot <- list(name = names(suit1), RPO_sp1 = sp1,
                        RPO_sp2 = sp2, RPO_combined = combined)
        names(res_tot)[2:3] <- c(spec1, spec2)
        if (plot) {
          pal <- c("white", "#7FCDBB", "#41B6C4", "#2C7FB8",
                   "#00012E")
          dp1 <- as.data.frame(res_tot$RPO_combined, xy = TRUE)
          p1 <- ggplot2::ggplot(co1xx) + ggplot2::geom_sf(colour = "transparent",
                                                          pch = 21, fill = "transparent", size = 5) +
            ggplot2::geom_tile(data = dp1, ggplot2::aes(x = dp1$x,
                                                        y = dp1$y, fill = dp1$RPO_combined)) + ggplot2::scale_fill_gradientn(colours = pal) +
            ggplot2::theme(panel.background = ggplot2::element_rect(fill = "gray90",
                                                                    colour = "black"), axis.title = ggplot2::element_blank(),
                           axis.text = ggplot2::element_text(size = 13),
                           legend.position = "none")
          dp2 <- as.data.frame(res_tot$RPO_combined_nok,
                               xy = TRUE)
          p2 <- ggplot2::ggplot(co1xx) + ggplot2::geom_sf(colour = "transparent",
                                                          pch = 21, fill = "transparent", size = 5) +
            ggplot2::geom_tile(data = dp2, ggplot2::aes(x = dp2$x,
                                                        y = dp2$y, fill = dp2$RPO_combined_nok)) +
            ggplot2::scale_fill_gradientn(colours = pal) +
            ggplot2::theme(panel.background = ggplot2::element_rect(fill = "gray90",
                                                                    colour = "black"), axis.title = ggplot2::element_blank(),
                           axis.text = ggplot2::element_text(size = 13),
                           legend.position = "none")
          pp <- cowplot::plot_grid(p1, p2, ncol = 2, align = "hv")
          pp
        }
        writeRaster(res_tot[[2]], paste(out.dir, "/", paste(spec1,
                                                            spec2, sep = "_"), "/", res_tot$name, "_", spec1,
                                        ".tif", sep = ""), overwrite = TRUE)
        writeRaster(res_tot[[3]], paste(out.dir, "/", paste(spec1,
                                                            spec2, sep = "_"), "/", res_tot$name, "_", spec2,
                                        ".tif", sep = ""), overwrite = TRUE)
        writeRaster(res_tot[[4]], paste(out.dir, "/", paste(spec1,
                                                            spec2, sep = "_"), "/", res_tot$name, "_combined.tif",
                                        sep = ""), overwrite = TRUE)
      }else res_tot <- NULL
    } else res_tot <- NULL
    return(res_tot)
  }
  outcome <- lapply(1:nlyr(pred[[1]]), function(layer) {
    RRphylogeo_inner(spec1, spec2, pred = lapply(pred, "[[",
                                                 layer), aggr = aggr, occs = occs, time_col = time_col,
                     th = th, kde_inversion = kde_inversion, resistance_map = resistance_map,
                     mask_for_pred = mask_for_pred, clust = clust, standardize = standardize,
                     plot = plot)
  })
  return(outcome)
}
