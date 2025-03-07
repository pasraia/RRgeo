#' @importFrom terra focalMat distance as.polygons
#' @importFrom sf st_convex_hull st_union st_bbox
occ.desaggregation.RASTER<-function (df, colxy, rast, plot = TRUE) {
  df_ini <- df
  obj <- cellFromXY(rast, df[, colxy])
  if (any(is.na(obj))) {
    stop("no NA admitted in species occurrences!")
  }
  obj1 <- split(obj, obj)
  l <- sapply(obj1, length)
  if (max(l) > 1) {
    l1 <- l[l > 1]
    for (j in 1:length(l1)) {
      w <- which(obj %in% as.numeric(names(l1[j])))
      found <- df[w, ]
      s <- sample(w, 1)
      df[w[w != s], colxy] <- NA
    }
    df_final <- df[!is.na(df[, colxy[1]]), ]
    if (plot == TRUE) {
      dev.new()
      plot(df_ini[, colxy], main = "distribution of occurences",
           sub = paste("# initial (black):", nrow(df_ini),
                       " | # kept (red): ", nrow(df_final)), pch = 19,
           col = "black", cex = 0.5)
      points(df_final[, colxy], pch = 19, col = "red",
             cex = 0.2)
    }
    return(df_final)
  }
  if (max(l) == 1)  return(df)
}

fix.coastal.points<-function (data,r, ncell, occ.desaggregation){
  init <- nrow(data)
  st_crs(data)->crs_data
  cbind(st_coordinates(data),st_drop_geometry(data))->data
  xy.cols=1:2
  cc <- extract(r,data[,xy.cols],ID = FALSE,xy=TRUE,na.rm=FALSE)
  vv <- which(is.na(cc[, 1]))
  if (length(vv) != 0) {
    m <- focalMat(r, type = "rectangle", d = res(r)[1] * ncell)
    m[] <- 1
    m[ncell + 1, ncell + 1] <- 0
    tt <- r
    test <- pblapply(vv, function(k) {
      focal <- rbind(data[k, xy.cols])
      ad <- adjacent(r, cellFromXY(r, focal), directions = m,
                     pairs = FALSE)
      valori <- tt[ad[1, ]]
      new.coords <- rbind(xyFromCell(r, ad[1, ])[!is.na(valori[,
                                                               1]), ])
      focal <- vect(as.matrix(focal), crs = crs(r))
      new.coords <- vect(new.coords, crs = crs(r))
      if (nrow(new.coords) != 0) {
        new1 <- distance(new.coords, focal)
        minimo <- which.min(new1)
        ret <- new.coords[minimo[sample(1:length(minimo),
                                        1)], ]
        return(crds(ret))
      }else {
        return(c(NA, NA))
      }
    })
    test <- do.call(rbind, test)
    final1 <- data
    final1[vv, xy.cols] <- test
    final1 <- final1[complete.cases(final1[, xy.cols]), ]
    if (occ.desaggregation & nrow(final1) > 0) {
      final2 <- occ.desaggregation.RASTER(final1, xy.cols,r, F)
      final2 <-st_as_sf(final2,coords=xy.cols,crs=crs_data)

    }else {
      final2 <- final1
      final2 <-st_as_sf(final2,coords=xy.cols,crs=crs_data)
    }
    if (occ.desaggregation == FALSE)
      final2 <- final1
    cat("\n", "initial sample size: ", init, sep = "")
    cat("\n", "shifted points: ", sum(!is.na(test[, 1])),
        sep = "")
    cat("\n", "removed points: ", sum(is.na(test[, 1])),
        sep = "")
    if (occ.desaggregation)
      cat("\n", "removed duplicates: ", nrow(final1) -
            nrow(final2), sep = "")
    cat("\n", "final sample size: ", nrow(final2), sep = "")
    cat("\n")
    final2 <-st_as_sf(final2,coords=xy.cols,crs=crs_data)
    return(final2)
  }
  if (length(vv) == 0) {
    if (occ.desaggregation) {
      final2 <- occ.desaggregation.RASTER(data, xy.cols,r, FALSE)
      final2 <-st_as_sf(final2,coords=xy.cols,crs=crs_data)
      return(final2)
    }else return(data)
  }
}

density_background<-function(pres.locs, MASK, rm.pres=TRUE) {
  test <- as.polygons(MASK)
  test1 <- st_as_sf(test)
  test2 <- spatstat.geom::as.owin(st_bbox(test1))
  ppp <- spatstat.geom::ppp(st_coordinates(pres.locs)[, 1],
                            st_coordinates(pres.locs)[, 2],
                            window = test2,
                            checkdup=FALSE)
  dens <- spatstat.explore::density.ppp(ppp, dimyx = c(nrow(MASK),
                                                       ncol(MASK)))
  mm1 <- MASK
  mm1[!is.na(mm1)] <- 1
  dens.ras <- resample(rast(dens), mm1)
  dens.ras <- dens.ras * mm1
  dens.ras[] <- scales::rescale(values(dens.ras))
  if (rm.pres) dens.ras <- mask(dens.ras, pres.locs, inverse = TRUE)
  return(dens.ras)
}

utils::globalVariables(c("intcal20","shcal20", "marine20","model_outputs"))
