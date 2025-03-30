#'@title Project the ENFA and ENphylo models into new geographical space and
#'  time interval
#'@description The function projects species marginality and specialization
#'  factors in different geographical areas and timescales. The function is able
#'  to convert marginality and specialization factors in habitat suitability
#'  values by using the Mahalanobis distances method.
#'@usage ENphylo_prediction(object, newdata,
#'  convert.to.suitability=FALSE,output.dir,proj_name="outputs")
#'@param object a \code{list} of ENFA and ENphylo models.  Each element of the
#'  list must be named using the names of the modelled species.
#'@param newdata a \code{SpatRaster} object including explanatory variables onto
#'  which ENFA or ENphylo models are to be projected. The list of variables must
#'  match the list used to model the species.
#'@param convert.to.suitability logical. If \code{TRUE},
#'  \code{ENphylo_prediction} projects ENFA or ENphylo model predictions in
#'  different geographical areas and timescales.
#'@param output.dir the file path wherein \code{ENphylo_prediction} creates an
#'  "ENphylo_prediction" folder to store prediction outputs for each species.
#'@param proj_name name of the subfolder created within the individual species
#'  folders to contain the \code{ENphylo_prediction} outputs.
#'@author Alessandro Mondanaro, Mirko Di Febbraro, Silvia Castiglione, Carmela
#'  Serio, Marina Melchionna, Pasquale Raia
#'@details If \code{convert.to.suitability} is set as \code{TRUE},
#'  \code{ENphylo_prediction} uses the function
#'  \code{\link[adehabitatHS]{mahasuhab}} from the \pkg{adehabitatHS} R package
#'  (\cite{Calenge, 2006}) to compute the habitat suitability map of the species
#'  over a given area. The conversion of Mahalanobis distances into
#'  probabilities follows the chi-squared distribution. Specifically, we set the
#'  degree of freedom equal to \emph{n} rather than \emph{n-1} following
#'  \cite{Etherington (2019)}. To convert habitat suitability values into binary
#'  presence/absence values, \code{ENphylo_prediction} relies on three different
#'  thresholding methods available in the function
#'  \code{\link[PresenceAbsence]{optimal.thresholds}} (\cite{Freeman &
#'  Moisen, 2008}).
#'@seealso  \code{vignette("ENphylo", package = "RRgeo")}
#'@return The function stores all the results in a number of nested subfolders
#'  all contained in the "ENphylo_prediction" folder created in
#'  \code{output.dir}. This contains a subfolder for each individual species in
#'  \code{object}, in which a subfolder named according to \code{proj_name}
#'  contains all the outputs. Specifically, the function saves the predictions
#'  for marginality and specificity (more than one depending on the number of
#'  significant axes selected by \code{\link{ENphylo_modeling}}) in the new
#'  geographical areas along with the suitability and binarized maps. The latter
#'  are calculated by using the three different predefined thresholds:
#'  MaxSensSpec (i.e. maximize TSS), SensSpec (i.e. equalize sensitivity and
#'  specificity) and 10th percentile of predicted probability.
#'@seealso \link{getENphylo_results}; \href{../doc/ENphylo.html}{\code{ENphylo} vignette}
#'@importFrom biomod2 bm_BinaryTransformation
#'@importFrom stats quantile
#'@importFrom terra app rast nlyr
#'@export
#'@references Calenge, C. (2006) The package adehabitat for the R software: a
#'  tool for the analysis of space and habitat use by animals. \emph{Ecological
#'  Modelling}, 197, 516-519.
#'@references Etherington, T. R. (2019). Mahalanobis distances and ecological
#'  niche modelling: correcting a chi-squared probability error. \emph{PeerJ},
#'  7, e6678.
#'@references Freeman, E. A. & Moisen, G. (2008). PresenceAbsence: An R Package
#'  for Presence-Absence Model Analysis. \emph{Journal of Statistical Software},
#'  23(11):1-31.
#'@references Mondanaro, A., Di Febbraro, M., Castiglione, S., Melchionna, M.,
#'  Serio, C., Girardi, G., Blefiore, A.M., & Raia, P. (2023). ENphylo: A new
#'  method to model the distribution of extremely rare species. \emph{Methods in
#'  Ecology and Evolution}, 14: 911-922. doi:10.1111/2041-210X.14066
#'@examples
#' \donttest{
#' library(ape)
#' library(terra)
#' library(sf)
#' library(RRgeo)
#'
#' newwd<-tempdir()
#' # newwd<-"YOUR_DIRECTORY"
#'
#' latesturl<-RRgeo:::get_latest_version("12734585")
#' curl::curl_download(url = paste0(latesturl,"/files/dat.Rda?download=1"),
#'                     destfile = file.path(newwd,"dat.Rda"), quiet = FALSE)
#' load(file.path(newwd,"dat.Rda"))
#' read.tree(system.file("exdata/Eucopdata_tree.txt", package="RRgeo"))->tree
#' tree$tip.label<-gsub("_"," ",tree$tip.label)
#' curl::curl_download(paste0(latesturl,"/files/X35kya.tif?download=1"),
#'                     destfile = file.path(newwd,"X35kya.tif"), quiet = FALSE)
#' rast(file.path(newwd,"X35kya.tif"))->map35
#' project(map35,st_crs(dat[[1]])$proj4string,res = 50000)->map
#'
#' ENphylo_modeling(input_data=dat[c(1,11)],
#'                  tree=tree,
#'                  input_mask=map[[1]],
#'                  obs_col="OBS",
#'                  time_col="age",
#'                  min_occ_enfa=15,
#'                  boot_test_perc=20,
#'                  boot_reps=10,
#'                  swap.args=list(nsim=5,si=0.2,si2=0.2),
#'                  eval.args=list(eval_metric_for_imputation="AUC",
#'                                 eval_threshold=0.7,
#'                                 output_options="best"),
#'                  clust=NULL,
#'                  output.dir=newwd)
#'
#'
#' getENphylo_results(input.dir =newwd,
#'                    mods="all",
#'                    species_name=names(dat)[c(1,11)])->mod
#'
#'
#' library(rnaturalearth)
#' ne_countries(returnclass = "sf")->globalmap
#' subset(globalmap,continent=="North America")->ame_map
#'
#' map35[[c("bio1","bio4","bio11","bio19")]]->newmap
#' crop(newmap,ext(ame_map))->newmap
#' project(newmap,st_crs(dat[[1]])$proj4string,res = 50000)->newmap
#'
#' ENphylo_prediction(object = mod,
#'                    newdata = newmap,
#'                    convert.to.suitability = TRUE,
#'                    output.dir=newwd,
#'                    proj_name="proj_example")
#'}

ENphylo_prediction<-function (object,
                              newdata,
                              convert.to.suitability = FALSE,
                              output.dir,
                              proj_name = "outputs")
{
  if (!(extends(class(newdata), "SpatRaster") | is.matrix(newdata) |
        is.data.frame(newdata)))
    stop("Please, provide newdata as a SpatRaster object or a data.frame")
  x <- newdata
  if (is.null(names(object))) {
    if (object[[1]]$call == "enfa") {
      U <- object[[1]]$co
      f1 <- function(y) y %*% U
      if (extends(class(x), "SpatRaster")) {
        if (!(all(names(x) %in% rownames(U)) & all(rownames(U) %in%
                                                   names(x))))
          stop("Variable names in newdata must match with those used to calibrate models") else {
            x <- x[[match(rownames(U), names(x))]]
          }
        ras <- app(x, fun = f1)
        names(ras) <- names(object[[1]]$sf)
      }else {
        if (!(all(colnames(x) %in% rownames(U)) & all(rownames(U) %in%
                                                      colnames(x))))
          stop("Variable names in newdata must match with those used to calibrate models")else {
            x <- x[, match(rownames(U), colnames(x))]
          }
        ras <- t(apply(x, 1, f1))
        colnames(ras) <- names(object[[1]]$sf)
        ras <- as.data.frame(ras)
      }
    }
    if (object[[1]]$call == "calibrated_enfa") {
      U <- object[[1]]$full_model$co
      f1 <- function(y) y %*% U
      if (extends(class(x), "SpatRaster")) {
        if (!(all(names(x) %in% rownames(U)) & all(rownames(U) %in%
                                                   names(x))))
          stop("Variable names in newdata must match with those used to calibrate models")else {
            x <- x[[match(rownames(U), names(x))]]
          }
        ras <- app(x, fun = f1)
        names(ras) <- names(object[[1]]$full_model$sf)
      }else {
        if (!(all(colnames(x) %in% rownames(U)) & all(rownames(U) %in%
                                                      colnames(x))))
          stop("Variable names in newdata must match with those used to calibrate models")else {
            x <- x[, match(rownames(U), colnames(x))]
          }
        ras <- t(apply(x, 1, f1))
        colnames(ras) <- names(object[[1]]$full_model$sf)
        ras <- as.data.frame(ras)
      }
    }
    ras <- list(call = "enfa_prediction", enfa_prediction = ras)
  }else {
    if(is.null(output.dir))  stop('argument "output.dir" is missing, with no default')
    message(paste("\n", "PREDICTING ENFA/IMPUTATION", "\n"))
    ras <- pblapply(object, function(ob) {
      if (ob$call == "calibrated_enfa") {
        U <- ob$calibrated_model$full_model$co
        f1 <- function(y) y %*% U
        if (!(all(names(x) %in% rownames(U)) & all(rownames(U) %in%
                                                   names(x))))
          stop("Variable names in newdata must match with those used to calibrate models") else {
            x <- x[[match(rownames(U), names(x))]]
          }
        ras <- app(x, fun = f1)
        names(ras) <- names(ob$calibrated_model$full_model$sf)
        ras <- ras[[1:ob$calibrated_model$full_model$significant_axes]]
        ras
        ras <- list(call = "enfa_prediction", enfa_prediction = ras)
      }
      if (ob$call == "calibrated_imputed") {
        ras <- lapply(ob$calibrated_model$co, function(co) {
          U <- co
          f1 <- function(y) y %*% U
          if (!(all(names(x) %in% rownames(U))&
                all(rownames(U)%in%names(x))))
            stop("Variable names in newdata must match with those used to calibrate models") else
              x <- x[[match(rownames(U), names(x))]]

          ras <- app(x, fun = f1)
          names(ras) <- colnames(U)
          ras
        })
        ras <- lapply(1:nlyr(ras[[1]]), function(i) rast(lapply(ras,
                                                                "[[", i)))
        names(ras)[1] <- "Marg"
        names(ras)[-1] <- paste("Spec", 1:(length(ras) -
                                             1), sep = "")

        if (ob$calibrated_model$output_options[1] ==
            "weighted.mean" & ob$calibrated_model$output_options[2] ==
            "OMR")
          warning("The weighted mean of model predictions using OMR provides misleading results")
        if (ob$calibrated_model$output_options[1] ==
            "weighted.mean") {
          ras <- c(Reduce("+", Map("*", lapply(1:nlyr(ras[[1]]),
                                               function(i) rast(lapply(ras, "[[", i))),
                                   ob$calibrated_model$evaluation[, ob$calibrated_model$output_options[2]]))/nlyr(ras[[1]]))

          lapply(1:nlyr(ras),function(ll){
            ras[[ll]]->cc
            names(cc)<-"swap_weighted"
            cc
          })->cc
          names(cc)<-names(ras)
          cc->ras
        }

        if (ob$calibrated_model$output_options[1] !=
            "weighted.mean") {

          ras <- lapply(ras, function(o) {

            if (any(grepl("evaluation", names(ob$calibrated_model)))) {
              names(o) <- rownames(ob$calibrated_model$evaluation)
            }else {
              names(o) <- paste("swap", 1:nlyr(o), sep = "_")
            }
            o
          })
        }
        ras <- list(call = "imputed_prediction", imputed_prediction = ras)
      }
      return(ras)
    })
    #setwd(output.dir)
    # suppressWarnings({dir.create("ENphylo_prediction")})
    ENout.dir<-file.path(output.dir,"ENphylo_prediction")
    dir.create(ENout.dir,showWarnings = FALSE)
    lapply(1:length(ras), function(xx) {
      dir.create(file.path(ENout.dir, names(ras)[xx],
                           proj_name), recursive = TRUE)
      if (ras[[xx]]$call == "imputed_prediction") {
        # dir.create(paste0("ENphylo_prediction/", names(ras)[xx],
        #                   "/", proj_name), recursive = TRUE)
        # lapply(1:length(ras[[xx]]$imputed_prediction),function(jj){
        #   writeRaster(ras[[xx]]$imputed_prediction[[jj]],
        #               file.path("ENphylo_prediction", names(ras)[xx],
        #                         proj_name,paste0(names(ras[[xx]]$imputed_prediction)[jj],
        #                                          ".tif")), overwrite = TRUE)
        # })


        lapply(1:length(ras[[xx]]$imputed_prediction),function(jj){
          writeRaster(ras[[xx]]$imputed_prediction[[jj]],
                      file.path(ENout.dir, names(ras)[xx],
                                proj_name,paste0(names(ras[[xx]]$imputed_prediction)[jj],
                                                 ".tif")), overwrite = TRUE)
        })
      }
      if (ras[[xx]]$call == "enfa_prediction") {
        # dir.create(paste0("ENphylo_prediction/", names(ras)[xx],
        #                   "/", proj_name), recursive = TRUE)
        # lapply(1:nlyr(ras[[xx]]$enfa_prediction), function(jj) {
        #   writeRaster(ras[[xx]]$enfa_prediction[[jj]],
        #               file.path("ENphylo_prediction", names(ras)[xx],
        #                         proj_name, paste0(names(ras[[xx]]$enfa_prediction[[jj]]),
        #                                           ".tif")), overwrite = TRUE)
        # })
        # dir.create(file.path(ENout.dir, names(ras)[xx],
        #                      proj_name), recursive = TRUE)
        lapply(1:nlyr(ras[[xx]]$enfa_prediction), function(jj) {
          writeRaster(ras[[xx]]$enfa_prediction[[jj]],
                      file.path(ENout.dir, names(ras)[xx],
                                proj_name, paste0(names(ras[[xx]]$enfa_prediction[[jj]]),
                                                  ".tif")), overwrite = TRUE)
        })
      }
    })
  }
  if (convert.to.suitability) {
    message(paste("\n", "CONVERTING PREDICTED VALUES TO SUITABILITY",
              "\n"))
    reference <- as.data.frame(x)
    ras_suitability <- pblapply(names(object), function(sp) {
      message(paste("\n", sp, "\n"))
      mydata <- object[[sp]]$formatted_data
      obs_col <- mydata$obs_col
      time_col <- mydata$time_col
      geoID_col <- mydata$geoID_col
      reference$type <- factor("reference", levels = c("obs_background",
                                                       "reference"))
      reference[, obs_col] <- 0
      if (nrow(mydata$input_back) > 10000) {
        independent_data_for_pred <- rbind(mydata$input_ones,
                                           mydata$input_back[sample(nrow(mydata$input_back),
                                                                    10000), ])
      }else independent_data_for_pred <- rbind(mydata$input_ones,mydata$input_back)

      obs_background <- independent_data_for_pred
      obs_background <- obs_background[, !grepl(geoID_col,
                                                colnames(obs_background))]
      obs_background$type <- factor("obs_background", levels = c("obs_background",
                                                                 "reference"))
      obs_background <- obs_background[, match(colnames(reference),
                                               colnames(obs_background))]
      ras_for_proj <- rbind(reference, obs_background)
      ras_obs <- ras_for_proj[, obs_col, drop = FALSE]
      ras_type <- ras_for_proj[, "type", drop = FALSE]
      ras_for_proj <- ras_for_proj[, !grepl(paste(c("type",
                                                    obs_col), collapse = "|"), colnames(ras_for_proj))]
      if (object[[sp]]$call == "calibrated_enfa") {
        U <- object[[sp]]$calibrated_model$full_model$co
        f1 <- function(y) y %*% U
        ma1 <- t(apply(as.matrix(ras_for_proj), 1, f1))
        colnames(ma1) <- colnames(U)
        ma1 <- as.data.frame(scale(ma1[, 1:object[[sp]]$calibrated_model$full_model$significant_axes]))
        ma1 <- list(ma1)
      }
      if (object[[sp]]$call == "calibrated_imputed") {
        ma1 <- lapply(1:length(object[[sp]]$calibrated_model$co),
                      function(kk) {
                        U <- object[[sp]]$calibrated_model$co[[kk]]
                        f1 <- function(y) y %*% U
                        ma1 <- t(apply(as.matrix(ras_for_proj), 1,
                                       f1))
                        colnames(ma1) <- colnames(U)
                        ma1 <- as.data.frame(scale(ma1))
                      })
        if (object[[sp]]$calibrated_model$output_options[1] ==
            "weighted.mean") {
          ma1 <- list(Reduce("+", Map("*", ma1, object[[sp]]$calibrated_model$evaluation[,
                                                                                         object[[sp]]$calibrated_model$output_options[2]]))/length(ma1))
        }
        else ma1 <- ma1
      }
      suitability_final <- lapply(1:length(ma1), function(kk) {
        repeat {
          hsm1 <- suppressWarnings(try(mahasuhab.custom(x = ma1[[kk]],
                                                        pts = ma1[[kk]][which(ras_obs == 1), ]),
                                       silent = TRUE))
          if (inherits(hsm1, "try-error")) {
            ma1[[kk]] <- ma1[[kk]][, -ncol(ma1[[kk]]),
                                   drop = FALSE]
            warning(paste("Dimensionality was reduced to",
                          ncol(ma1[[kk]]), "axes"))
          }
          else break
        }
        ma1[[kk]]$Suitability <- hsm1$MD
        reference <- ma1[[kk]][which(ras_type == "reference"),
        ]
        pred_for_th <- ma1[[kk]][which(ras_type == "obs_background"),
        ]
        data_for_th <- data.frame(ID = 1:nrow(pred_for_th),
                                  observed = obs_background[, obs_col], pred = pred_for_th$Suitability)
        suppressWarnings(TH <- optimal.thresholds(data_for_th)[2:3,
                                                               2])
        TH <- c(TH, quantile(data_for_th[which(data_for_th$observed ==
                                                 1), ]$pred, 0.1))
        hsm1 <- x[[1]]
        hsm1[!is.na(hsm1)] <- reference$Suitability
        names(hsm1) <- "Suitability"
        if (object[[sp]]$call == "calibrated_imputed") {
          if (any(grepl("evaluation", names(object[[sp]]$calibrated_model)))) {
            if (object[[sp]]$calibrated_model$output_options[1] ==
                "best") {
              names(hsm1) <- paste("Suitability", rownames(object[[sp]]$calibrated_model$evaluation)[kk],
                                   sep = "_")
            }
            if (object[[sp]]$calibrated_model$output_options[1] ==
                "weighted.mean") {
              names(hsm1) <- paste(names(hsm1), "swap_weighted",
                                   sep = "_")
            }
            if (object[[sp]]$calibrated_model$output_options[1] ==
                "full") {
              names(hsm1) <- paste("Suitability", "swap",
                                   kk, sep = "_")
            }
          }
          else names(hsm1) <- paste("Suitability", "swap",
                                    kk, sep = "_")
        }
        return(list(TH, hsm1))
      })
      if (object[[sp]]$call == "calibrated_enfa") {
        thh <- lapply(suitability_final, function(xx) {
          xx <- rast(sapply(xx[[1]], function(yy) bm_BinaryTransformation(xx[[2]],
                                                                          yy)))
          names(xx) <- c("SensSpec", "MaxSensSpec", "TenPerc")
          xx
        })
        names(thh[[1]]) <- paste("Binary", names(thh[[1]]),
                                 sep = "_")
        ras[[sp]]$enfa_prediction <- c(ras[[sp]]$enfa_prediction,
                                       suitability_final[[1]][[2]], thh[[1]])
        ras_final <- list(call = "enfa_prediction", enfa_prediction = ras[[sp]]$enfa_prediction)
      }
      if (object[[sp]]$call == "calibrated_imputed") {
        thh <- lapply(suitability_final, function(xx) {
          xx <- rast(sapply(xx[[1]], function(yy) bm_BinaryTransformation(xx[[2]],
                                                                          yy)))
          names(xx) <- c("SensSpec", "MaxSensSpec", "TenPerc")
          xx
        })
        thh <- lapply(1:3, function(i) rast(lapply(thh,
                                                   "[[", i)))
        if (any(grepl("evaluation", names(object[[sp]]$calibrated_model)))) {
          if (object[[sp]]$calibrated_model$output_options[1] ==
              "best") {
            thh <- lapply(thh, function(jj) {
              names(jj) <- paste("Binary", strsplit(names(jj[[1]]),
                                                    "[.]")[[1]][1], rownames(object[[sp]]$calibrated_model$evaluation),
                                 sep = "_")
              jj
            })
          }
          if (object[[sp]]$calibrated_model$output_options[1] ==
              "weighted.mean") {
            thh <- lapply(thh, function(jj) {
              names(jj) <- paste("Binary", strsplit(names(jj[[1]]),
                                                    "[.]")[[1]][1], "swap_weighted", sep = "_")
              jj
            })
          }
          if (object[[sp]]$calibrated_model$output_options[1] ==
              "full") {
            thh <- lapply(thh, function(jj) {
              names(jj) <- paste(strsplit(names(jj[[1]]),
                                          "[.]")[[1]][1], "swap", 1:nlyr(jj), sep = "_")
              jj
            })
          }
        }
        else {
          thh <- lapply(thh, function(jj) {
            names(jj) <- paste(strsplit(names(jj[[1]]),
                                        "[.]")[[1]][1], "swap", 1:nlyr(jj), sep = "_")
            jj
          })
        }
        if (object[[sp]]$calibrated_model$output_options[1] ==
            "full") {
          ras[[sp]]$imputed_prediction <- c(ras[[sp]]$imputed_prediction,
                                            Suitability = rast(lapply(suitability_final,
                                                                      "[[", 2)), Binary_SensSpec = thh[[1]],
                                            Binary_MaxSensSpec = thh[[2]], Binary_TenPerc = thh[[3]])
        }
        else {
          ras[[sp]]$imputed_prediction <- c(ras[[sp]]$imputed_prediction,
                                            Suitability = rast(lapply(suitability_final, "[[", 2)),
                                            Binary_SensSpec = thh[[1]],
                                            Binary_MaxSensSpec = thh[[2]],
                                            Binary_TenPerc = thh[[3]])
        }
        ras_final <- list(call = "imputed_prediction",
                          imputed_prediction = ras[[sp]]$imputed_prediction)
      }
      return(ras_final)
    })
    names(ras_suitability) <- names(ras)
    ras <- ras_suitability
    lapply(1:length(ras), function(xx) {
      if (ras[[xx]]$call == "imputed_prediction") {
        good <- grep("Suitability|Binary", names(ras[[xx]][[2]]))
        lapply(good, function(jj) {
          writeRaster(ras[[xx]][[2]][[jj]], file.path(ENout.dir,
                                                      names(ras)[xx],
                                                      proj_name,
                                                      paste0(names(ras[[xx]][[2]])[jj],
                                                             ".tif")), overwrite = TRUE)
        })
      }
      if (ras[[xx]]$call == "enfa_prediction") {
        good <- grep("Suitability|Binary", names(ras[[xx]][[2]]))
        lapply(good, function(jj) {
          writeRaster(ras[[xx]][[2]][[jj]], file.path(ENout.dir,
                                                      names(ras)[xx], proj_name, paste0(names(ras[[xx]][[2]])[[jj]],
                                                                                        ".tif")), overwrite = TRUE)
        })
      }
    })
  }
  return(ras)
}


