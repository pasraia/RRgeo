#' @importFrom terra ncol nrow crs global ncell xyFromCell
#' @importFrom ks Hpi Hscv.diag kde
#' @importFrom sf st_geometry_type st_point
#' @importFrom stats runif var
sf.kde.mod<-function (x, y = NULL, bw = NULL, ref = NULL, res = NULL, standardize = FALSE,
                      scale.factor = NULL, mask = FALSE)
{
  if (missing(x))
    stop("x argument must be provided")
  ref.flag = inherits(ref, "SpatRaster")
  if (!inherits(x, c("sf", "sfc")))
    stop(deparse(substitute(x)), " must be a sf, or sfc object")
  if (unique(as.character(st_geometry_type(x))) != "POINT")
    stop(deparse(substitute(x)), " must be single-part POINT geometry")
  if (is.null(ref)) {
    if (!is.null(res)) {
      ref <- rast(ext(x), resolution = res)
    }else {
      ref <- rast(ext(x))
      message("defaulting to ", res(ref)[1], "x", res(ref)[2],
              " cell resolution")
    }
  }
  if (inherits(ref, "numeric")) {
    if (length(ref) != 4)
      stop("Need xmin, xmax, ymin, ymax bounding coordinates")
    if (!is.null(res)) {
      ref <- rast(ext(ref), resolution = res)
    }else {
      ref <- rast(ext(ref))
      message("defaulting to ", res(ref)[1], "x", res(ref)[2],
              " cell resolution")
    }
  }else {
    if (!inherits(ref, "SpatRaster"))
      stop(deparse(substitute(ref)), " must be a terra SpatRast object")
  }
  n <- c(nrow(ref), ncol(ref))

  eps<-1e-3
  st_coordinates(x)->coords
  for (i in 1:2) {
    if (var(coords[,i])==0){
      coords[,i] <- coords[,i]+runif(nrow(coords),-eps,eps)
    }
  }
  new_geom<-st_sfc(
    lapply(seq_len(nrow(coords)), function(j) st_point(coords[j,])),
    crs=st_crs(x))
  sf::st_geometry(x)<-new_geom

  if (is.null(bw)) {
    if (inherits(try(suppressWarnings(Hpi(st_coordinates(x)[,
                                                            1:2])), silent = TRUE), "try-error")) {
      bw <- min(res(ref))
      bw <- diag(bw^2, 2) + diag(eps, 2)
      message("Using automatic bandwidth: ")
    }else {
      bw = suppressWarnings(Hpi(st_coordinates(x)[,
                                                  1:2]))
      res_min <- min(res(ref))
      if (any(diag(bw)<res_min^2)){
        bw<-diag(pmax(diag(bw),res_min^2))
      }

      message("Using automatic bandwidth: ")
    }
  }else {
    message("Using specified bandwidth: ")
  }
  if (!is.null(y)) {
    message("\n", "calculating weighted kde", "\n")
    kde.est <- suppressWarnings(rast(matrix(kde(st_coordinates(x)[,
                                                                  1:2], h = bw, eval.points = xyFromCell(ref, 1:ncell(ref)),
                                                gridsize = n, w = y, density = TRUE)$estimate, nrow = n[1],
                                            ncol = n[2], byrow = TRUE), extent = ext(ref)))
  }else {
    message("\n", "calculating unweighted kde", "\n")
    kde.est <- suppressWarnings(rast(matrix(kde(st_coordinates(x)[,
                                                                  1:2], H = bw, eval.points = xyFromCell(ref, 1:ncell(ref)),
                                                gridsize = n, density = TRUE)$estimate, nrow = n[1],
                                            ncol = n[2], byrow = TRUE), extent = ext(ref)))
  }
  if (!is.null(scale.factor))
    kde.est <- kde.est * scale.factor
  if (standardize == TRUE) {
    kde.est <- kde.est/global(kde.est, "max", na.rm = TRUE)[,
                                                            1]
  }
  if (mask) {
    if (!ref.flag) {
      message("Since a raster was not used as ref, there is nothing to mask")
    }
    else {
      kde.est <- mask(kde.est, ref)
    }
  }
  terra::crs(kde.est) <- crs(x)
  return(kde.est)
}
