#'@title Radiocarbon Calibration of Occurrences
#'@description The function is meant to automatically apply the calibration
#'  process to conventional radiocarbon ages, relying on the package Bchron
#'  (Haslett & Parnell 2008). Specifically, the function internally selects the
#'  appropriate calibration curve based on the latitude associated with each
#'  occurrence and the nature of the sample (i.e marine or terrestrial samples).
#'@usage cal14C(dataset,age=NULL,uncertainty=NULL,latitude=NULL,domain=NULL,
#'  bounds=c(0.025,0.975), clust=0.5, save=TRUE, output.dir=".")
#'@param dataset a \code{data.frame} containing all the occurrences to be
#'  calibrated.
#'@param age character. Name of the column in \code{dataset} containing the
#'  conventional radiocarbon dates.
#'@param uncertainty character. Name of the column in \code{dataset} containing
#'  the uncertainty associated to conventional radiocarbon dates.
#'@param latitude character. Name of the column in \code{dataset} containing the
#'  latitude (in decimal degrees) of each occurrence.
#'@param domain character. Name of the column in \code{dataset} indicating which
#'  occurrences are marine radiocarbon samples. If \code{NULL}, all the
#'  occurrences are assumed as "terrestrial" radiocarbon samples.
#'@param bounds numeric. An upper and lower bound (in quantiles) to define the
#'  limits of the density probability created for each radiocarbon age (default:
#'  95\%).
#'@param clust numeric. The proportion of cores used to train \code{cal14C}. If
#'  \code{NULL}, parallel computing is disabled.
#'@param save,output.dir if \code{save = TRUE}, \code{cal14C} outputs are saved
#'  in \code{output.dir}.
#'@author Alessandro Mondanaro, Silvia Castiglione, Pasquale Raia
#'@details If \code{dataset} includes marine samples, the user should indicate
#'  it in the domain column by indicating "marine" for the corresponding
#'  occurrences. In this case, the function uses the marine20 curve to calibrate the
#'  related radiocarbon ages accordingly.
#'@return The initial \code{dataset} with additional
#'  columns providing detailed calibration information for each occurrence. The
#'  new columns indicate the calibration curve used for each occurrence
#'  ("curve"), the calibrated radiocarbon ages ("cal_age"), and the values
#'  corresponding to the specified confidence limits derived from the density
#'  estimate of the calibrated radiocarbon ages. If \code{save=TRUE},
#'  the dataframe is saved as xlsx file in \code{output.dir}.
#'@importFrom parallel clusterEvalQ clusterExport
#'@importFrom utils data
#'@export
#'@seealso \link[Bchron]{BchronCalibrate}
#'@examples
#' \donttest{
#'
#' library(RRgeo)
#'
#' ## Create an example dataset with 100 random radiocarbon ages and errors
#' set.seed(2025)
#' data.frame(age=round(runif(100,20000,50000),0),
#'            uncertain=round(runif(100,20,300),0),
#'            latitude=round(runif(100,-90,90),2))->data
#' data$domain<-"domain"
#' rep("marine",5)-> data[sample(nrow(data),5),"domain"]
#'
#' cal14C(dataset=data,
#'        age="age",
#'        uncertainty = "uncertain",
#'        latitude = "latitude",
#'        domain<-"domain",
#'        clust=NULL,
#'        save= FALSE)->res
#'
#'}

cal14C<-function (dataset,
                      age = NULL,
                      uncertainty = NULL,
                      latitude = NULL,
                      domain = NULL,
                      bounds = c(0.025, 0.975),
                      clust = 0.5,
                      save = TRUE,
                      output.dir = "."){
  if (!requireNamespace("Bchron", quietly = TRUE))
    stop("Package \"Bchron\" needed for this function to work, please install it.",
         call. = FALSE)
  if (save && (!requireNamespace("openxlsx", quietly = TRUE))) {
    stop("Package \"openxlsx\" needed for save=TRUE. Please install it.",
         call. = FALSE)
  }
  dat <- dataset
  dat$ID <- seq(1:length(dat[, 1]))
  if (!is.null(age))
    colnames(dat)[which(colnames(dat) == age)] <- "age"
  if (!is.null(uncertainty))
    colnames(dat)[which(colnames(dat) == uncertainty)] <- "error"
  if (!is.null(latitude))
    colnames(dat)[which(colnames(dat) == latitude)] <- "lat"
  if (!is.null(domain))
    colnames(dat)[which(colnames(dat) == domain)] <- "domain"
  dat$age <- round(dat$age)
  dat$error <- round(dat$error)
  dat$curve <- NA
  if (any(which(dat$lat >= 0))) {
    data("intcal20", envir = environment(), package = "Bchron")
    dat[which(dat$lat >= 0 & dat$age >= range(intcal20$V2)[1] &
                dat$age <= range(intcal20$V2)[2]), ]$curve <- "intcal20"
  }
  if (any(which(dat$lat < 0))) {
    data("shcal20", envir = environment(), package = "Bchron")
    dat[which(dat$lat < 0 & dat$age >= range(shcal20$V2)[1] &
                dat$age <= range(shcal20$V2)[2]), ]$curve <- "shcal20"
  }
  if (!is.null(dat$domain)) {
    dat[which(dat$domain == "marine"), ]$curve <- NA
    data("marine20", envir = environment(), package = "Bchron")
    dat[which(dat$domain == "marine" & dat$age >= range(marine20$V2)[1] &
                dat$age <= range(marine20$V2)[2]), ]$curve <- "marine20"
  }
  if (any(which(is.na(dat$curve)))) {
    add_dat <- dat[which(is.na(dat$curve)), ]
    add_dat$curve <- "not calibrated"
    add_dat$cal.age <- add_dat$CI_97.5 <- add_dat$CI_2.5 <- NA
    colnames(add_dat)[grepl("CI_2.5|CI_97.5", colnames(add_dat))] <- paste0("CI_",
                                                                            c(bounds[1] * 100, bounds[2] * 100))
    dat <- dat[-which(is.na(dat$curve)), ]
  }else add_dat <- NULL

  message("\n", "PERFORMING CALIBRATION", "\n")
  if (!is.null(clust)) {
    cl <- detectCores() * clust
    cl <- makeCluster(cl)
    clusterExport(cl, "dat", envir = environment())
    clusterEvalQ(cl, {
      library(Bchron)
    })
    cal <- pblapply(1:length(dat$curve), function(x) {
      Bchron::BchronCalibrate(dat$age[x], dat$error[x],
                              calCurves = dat$curve[x])
    }, cl = cl)
    closeAllConnections()
    stopCluster(cl)
    gc()
  }else {
    cal <- pblapply(1:length(dat$curve), function(x) {
      Bchron::BchronCalibrate(dat$age[x], dat$error[x],
                              calCurves = dat$curve[x])
    })
  }
  mean.cal <- lapply(1:length(cal), function(k) mean(cal[[k]]$Date1$ageGrid))
  c95 <- lapply(1:length(cal), function(k) quantile(cal[[k]]$Date1$ageGrid,
                                                    c(bounds[1], bounds[2])))
  calibrated <- data.frame(dat, cal.age = round(do.call(rbind,
                                                        mean.cal)), CI_2.5 = round(do.call(rbind, c95))[, 1],
                           CI_97.5 = round(do.call(rbind, c95)[, 2]))
  colnames(calibrated)[grepl("CI_2.5|CI_97.5", colnames(calibrated))] <- paste0("CI_",
                                                                                c(bounds[1] * 100, bounds[2] * 100))

  add_dat[match(colnames(calibrated),colnames(add_dat))]->add_dat
  dat.cal <- rbind(calibrated, add_dat)
  dat.cal <- dat.cal[order(dat.cal$ID), ]
  dat.cal <- dat.cal[, -which(colnames(dat.cal) == "ID")]
  if (!is.null(age))
    colnames(dat.cal)[which(colnames(dat.cal) == "age")] <- age
  if (!is.null(uncertainty))
    colnames(dat.cal)[which(colnames(dat.cal) == "error")] <- uncertainty
  if (!is.null(latitude))
    colnames(dat.cal)[which(colnames(dat.cal) == "lat")] <- latitude
  if (!is.null(domain))
    colnames(dat.cal)[which(colnames(dat.cal) == "domain")] <- domain
  if (save)
    openxlsx::write.xlsx(dat.cal, paste(output.dir, "dat.cal.xlsx",
                                        sep = "/"))
  gc()
  return(dat.cal)
}

