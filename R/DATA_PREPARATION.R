#' @importFrom sf st_crs st_coordinates st_drop_geometry
#' @importFrom grDevices dev.new
#' @importFrom graphics points
#' @importFrom terra extract crs cellFromXY adjacent
#' @importFrom stats complete.cases
DATA_PREPARATION<-function(species_input_data,
                           obs_col = "OBS",
                           input_mask,
                           time_col = NULL){
  if (!crs(input_mask, proj = TRUE) == st_crs(species_input_data)$proj4string |
      is.na(crs(input_mask, proj = TRUE)) | is.na(st_crs(species_input_data)$proj4string))
    stop("some missing or not-matching projections in input data or background area")
  occ.desaggregation.RASTER <- function(df, colxy, rast, plot = TRUE) {
    df_ini <- df
    cac <- cellFromXY(rast, df[, colxy])
    if (any(is.na(cac))) {
      stop("no NA admitted in species occurrences!")
    }
    cac1 <- split(cac, cac)
    l <- sapply(cac1, length)
    if (max(l) > 1) {
      l1 <- l[l > 1]
      df_ok <- df[which(cac %in% as.numeric(names(l[l ==
                                                        1]))), ]
      found <- lapply(1:length(l1), function(j) {
        w <- which(cac %in% as.numeric(names(l1[j])))
        s <- sample(w, 1)
        df[s, ]
      })
      df_final <- rbind(df_ok, do.call(rbind, found))
      if (plot) {
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
    if (max(l) == 1)
      return(df)
  }
  species_input_data <- cbind(st_coordinates(species_input_data),
                              as.data.frame(st_drop_geometry(species_input_data)))
  species_input_data0 <- subset(species_input_data, species_input_data[,
                                                                       obs_col] == 0)
  species_input_data <- subset(species_input_data, species_input_data[,
                                                                      obs_col] == 1)
  if (is.null(time_col)) {
    species_input_data2 <- occ.desaggregation.RASTER(species_input_data,
                                                     1:2, rast = input_mask, plot = FALSE)
    if (nrow(species_input_data2) < nrow(species_input_data))
      warning(paste(nrow(species_input_data) - nrow(species_input_data2),
                    "duplicated points into raster cells were found and removed!"))
    species_input_data <- species_input_data2
  }else {
    ss <- split(species_input_data, species_input_data[,
                                                       time_col])
    species_input_data2 <- lapply(ss, function(x) {
      occ.desaggregation.RASTER(x, 1:2, rast = input_mask,
                                plot = FALSE)
    })
    species_input_data2 <- do.call(rbind, species_input_data2)
    if (nrow(species_input_data2) < nrow(species_input_data))
      warning(paste(nrow(species_input_data) - nrow(species_input_data2),
                    "duplicated points into raster cells were found and removed!"))
    species_input_data <- species_input_data2
  }
  species_input_data <- rbind(species_input_data, species_input_data0)
  colnames(species_input_data)[1:2] <- c("coordsX", "coordsY")
  coords_cols <- grep("coords", colnames(species_input_data))
  ones <- subset(species_input_data, species_input_data[, obs_col] ==
                   1)
  back <- subset(species_input_data, species_input_data[, obs_col] ==
                   0)
  input_ones <- cbind.data.frame(ones[, -coords_cols], geoID = extract(input_mask,
                                                                       ones[, coords_cols], cells = TRUE)[, "cell"])
  input_back <- cbind.data.frame(back[, -coords_cols], geoID = extract(input_mask,
                                                                       back[, coords_cols], cells = TRUE)[, "cell"])

  if (nrow(input_ones) < 5 && !is.null(time_col)) {
    split(input_ones,input_ones[,time_col])->aa
    split(input_back,input_back[,time_col])->bb
    split(ones,ones[,time_col])->aa1
    split(back,back[,time_col])->bb1

    lapply(1:length(aa),function(x){

      all_ones <- cbind(aa1[[x]][, coords_cols], aa[[x]])
      all_back <- cbind(bb1[[x]][, coords_cols], bb[[x]])
      ad <- adjacent(input_mask, all_ones$geoID, pairs = TRUE,
                     directions = 16)
      ad <- setdiff(unique(c(ad[, 2])), unique(c(ad[, 1])))
      focus <- all_ones
      ad <- all_back[all_back$geoID %in% ad, ]
      ad[complete.cases(ad), , drop = FALSE]

      focus[,!grepl(paste(c(obs_col, time_col, "geoID",
                            "coordsX", "coordsY"), collapse = "|"),
                    colnames(focus))]->ff
      focus <- focus[!duplicated(ff), ]
      ad[,!grepl(paste(c(obs_col, time_col, "geoID",
                         "coordsX", "coordsY"), collapse = "|"),
                 colnames(ad))]->aa
      ad <- ad[!duplicated(aa), ]

      if (any(duplicated(rbind(ff,aa)))) {
        ad <- ad[-which(as.numeric(duplicated(rbind(ff, aa)))[-c(1:nrow(ff))] ==
                          1), ]
      }
      list(focus=focus,ad=ad)
    })->all_data

    do.call(rbind,lapply(all_data,"[[",1))->good
    do.call(rbind,lapply(all_data,"[[",2))->to_add

    if (nrow(good) + nrow(to_add) >= 5 & nrow(to_add) > 0) {
      lapply(1:length(all_data),function(bb){
        ff<-all_data[[bb]]$focus
        aa<-all_data[[bb]]$ad

        lapply(1:nrow(ff), function(j) {
          mat <- as.matrix(rbind(ff[j, ], aa))
          lapply(2:nrow(mat), function(jj) {
            a<-mat[1,!grepl(paste(c(obs_col, time_col, "geoID",
                                    "coordsX", "coordsY"), collapse = "|"), colnames(mat))]
            b<-mat[jj,!grepl(paste(c(obs_col, time_col, "geoID",
                                     "coordsX", "coordsY"), collapse = "|"), colnames(mat)) ]

            ppA <- ((a %*% b)/(RRphylo:::unitV(a) * RRphylo:::unitV(b)))
            RRphylo:::rad2deg(acos(ppA))->vv
            cbind(val=vv,data.frame(t(mat[jj,])))
          })->closer
          do.call(rbind,closer)
        })->closer
        closer <- do.call(cbind, closer)
      })->all

      do.call(rbind,all)->all
      good <- all[order(all$val),][1:(5 - nrow(input_ones)),]
      good[,grep("val",colnames(good),invert=TRUE)]->good
      good[,grepl(obs_col,colnames(good))]<-1
      input_ones<-rbind(input_ones,good[!grepl(paste(c("coordsX",
                                                       "coordsY"), collapse = "|"), colnames(good))])
      ones<-rbind(ones,good[!grepl("geoID", colnames(good))])
    }
    else warning("Too few presence points: the model might be unreliable")
  }

  if (nrow(input_ones) < 5 && is.null(time_col)) {
    all_ones <- cbind(ones[, coords_cols], input_ones)
    all_back <- cbind(back[, coords_cols], input_back)
    ad <- adjacent(input_mask, all_ones$geoID, pairs = TRUE,
                   directions = 16)
    ad <- setdiff(unique(c(ad[, 2])), unique(c(ad[, 1])))
    focus <- all_ones
    ad <- all_back[all_back$geoID %in% ad, ]
    ad[complete.cases(ad), , drop = FALSE]
    focus[,!grepl(paste(c(obs_col, "geoID",
                          "coordsX", "coordsY"), collapse = "|"),
                  colnames(focus))]->ff
    focus <- focus[!duplicated(ff), ]
    ad[,!grepl(paste(c(obs_col, "geoID",
                       "coordsX", "coordsY"), collapse = "|"),
               colnames(ad))]->aa
    ad <- ad[!duplicated(aa), ]

    if (any(duplicated(rbind(ff,aa)))) {
      ad <- ad[-which(as.numeric(duplicated(rbind(ff, aa)))[-c(1:nrow(ff))] ==
                        1), ]
    }

    if (nrow(focus) + nrow(ad) >= 5 & length(ad) > 0) {
      lapply(1:nrow(focus), function(j) {
        mat <- as.matrix(rbind(focus[j, ], ad))
        lapply(2:nrow(mat), function(jj) {
          a<-mat[1,!grepl(paste(c(obs_col, "geoID",
                                  "coordsX", "coordsY"), collapse = "|"), colnames(mat))]
          b<-mat[jj,!grepl(paste(c(obs_col, "geoID",
                                   "coordsX", "coordsY"), collapse = "|"), colnames(mat)) ]


          ppA <- ((a %*% b)/(RRphylo:::unitV(a) * RRphylo:::unitV(b)))
          RRphylo:::rad2deg(acos(ppA))->vv
          cbind(val=vv,data.frame(t(mat[jj,])))
        })->closer
        do.call(rbind,closer)
      })->all
      all <- do.call(rbind, all)
      good <- all[order(all$val),][1:(5 - nrow(input_ones)),]
      good[,grep("val",colnames(good),invert=TRUE)]->good
      good[,grepl(obs_col,colnames(good))]<-1
      input_ones<-rbind(input_ones,good[!grepl(paste(c("coordsX",
                                                       "coordsY"), collapse = "|"), colnames(good))])
      ones<-rbind(ones,good[!grepl("geoID", colnames(good))])
    }
    else warning("Too few presence points: the model might be unreliable")
  }

  if (!all(input_ones$geoID %in% input_back$geoID))
    stop("some points fall outside of the background area")
  ones_coords <- ones[, coords_cols]
  rownames(ones_coords) <- NULL
  OUTPUT <- list(call = "input_data", input_ones = input_ones,
                 input_back = input_back, obs_col = obs_col, geoID_col = "geoID",
                 time_col = time_col, ones_coords = ones_coords, study_area = input_mask)
  return(OUTPUT)
}
