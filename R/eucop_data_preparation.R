#'@title Import and preprocess mammal occurrence data
#'@description The function is meant to automatically import and preprocess
#'  fossil mammal occurrences and paleoclimatic/vegetational data available in
#'  EutherianCop dataset (Mondanaro et al., 2025). It also provides two distinct
#'  approaches, both implemented within a user-defined study area, for sampling
#'  a specified number of pseudoabsences or alternatively defining the
#'  background points. This flexibility enables users to assemble a list of
#'  \code{sf} objects that can be easily used to train ENFA, ENphylo or any
#'  other SDM algorithms of their choice.
#'@usage
#'eucop_data_preparation(input.dir,species_name,variables="all",which.vars=NULL,
#'calibration=FALSE,add.modern.occs=FALSE,
#'combine.ages=NULL,remove.duplicates=TRUE, bk_points=NULL,output.dir=".")
#'@param input.dir the file path wherein EutherianCop mammal occurrences and
#'  paleoclimatic data are to be stored.
#'@param species_name character. The name of the single (or multiple) species
#'  used by \code{eucop_data_preparation}.
#'@param variables character. The name of paleoclimatic simulations to be used.
#'  The viable options are "climveg", "bio", or "all".
#'@param which.vars character vector indicating the name of the variables to be
#'  downloaded. The list of accepted names can be found
#'  [here](https://www.nature.com/articles/s41597-024-04181-4/tables/1).
#'@param calibration logical. If \code{TRUE}, \code{eucop_data_preparation}
#'  performs the 14C calibration process to convert the conventional radiocarbon
#'  age estimates included in EutherianCop raw data file.
#'@param add.modern.occs logical. If \code{TRUE}, \code{eucop_data_preparation}
#'  adds the modern records (if present) related to species in
#'  \code{species_name}.
#'@param combine.ages one of \code{"mean"} or \code{"median"}. The method to be
#'  used to aggregate multiple ages for each site or layer within the site.
#'@param remove.duplicates logical. If \code{TRUE},
#'  \code{eucop_data_preparation} removes duplicated record for each grid cell
#'  within a given time bin.
#'@param bk_points a list including parameters to add background/pseudoabsence
#'  (i.e. absence) points (following the procedure described in Mondanaro et al.
#'  2024). The list includes:
#'  \itemize{\item buff: the proportional distance to set a buffer around the
#'  minimum convex polygon that encompasses all occurrences of the target species.
#'  \item bk_strategy: the strategy to add the absence points. It can be one of
#'  "background" or "pseudoabsence". \item bk_n: number of absence points.}
#'  If provided as an empty \code{list()}, the function automatically sets
#'  \code{buff = 0.1}, \code{bk_strategy="background"},\code{bk_n=10000}.
#'@param output.dir the file path wherein \code{eucop_data_preparation} stores
#'  the results.
#'@author Alessandro Mondanaro, Silvia Castiglione, Pasquale Raia
#'@details The variables argument allows the selection of climatic and
#'  environmental variables ("climveg"), bioclimatic variables ("bio"), or both
#'  sets of variables.
#'
#'@details Through the \code{bk_strategy} argument,
#'  \code{eucop_data_preparation} offers two different approaches to generate
#'  absence points. The definition of the study area is the same for both
#'  methods. Under \code{bk_strategy = "background"}, the \code{bk_n} argument
#'  defines the maximum number of background points sampled from the study area
#'  within each time bin. Under \code{bk_strategy = "pseudoabsence"}, the
#'  \code{bk_n} argument represents the maximum number of pseudoabsence points
#'  across all time bins. This flexibility allows users to accommodate the
#'  different requirements for training the traditional envelope models (i.e.
#'  ENFA, ENphylo) and the common correlative or machine learning models (i.e.
#'  generalized linear model, MaxEnt, Random Forest).
#'@details Additionally, if \code{bk_points} is not \code{NULL}, the ages of
#'  presences and pseudoabsences or background points are forced to 1 kyr
#'  resolution according to the temporal resolution of the
#'  paleoclimatic/vegetational or bioclimatic data.
#'@return \code{eucop_data_preparation} does not store any results in the global
#'  environment. Instead, a list of GeoPackage files, one per selected species,
#'  is saved in the directory specified by \code{output.dir}. The names of these
#'  files depend on the combination of arguments chosen by users: they include
#'  the suffix "cal/uncal" and "combined/multi" depending on whether calibration
#'  (\code{calibration}) and age aggregation (\code{combine.ages}) steps are
#'  performed. In any case, output files include information about ages, a
#'  column called "OBS" including species occurrence data in binary format,
#'  spatial geometry, and all the data information derived from EutherianCop
#'  dataset.
#'@importFrom gtools mixedorder
#'@importFrom sf sf_use_s2 st_distance st_buffer st_crop st_join st_write
#'  st_convex_hull st_transform
#'@importFrom terra mask as.points project
#'@seealso \href{../doc/Preparing-Data.html}{\code{eucop_data_preparation} vignette}
#'@export
#'@references Mondanaro, A., Di Febbraro, M., Castiglione, S., Belfiore, A. M.,
#'  Girardi, G., Melchionna, M., Serio, C., Esposito, A., & Raia, P. (2024).
#'  Modelling reveals the effect of climate and land use change on Madagascar’s
#'  chameleons fauna. \emph{Communications Biology}, 7: 889.
#'  doi:10.1038/s42003-024-06597-5.
#'@references Mondanaro, A., Girardi, G., Castiglione, S., Timmermann, A.,
#'  Zeller, E., Venugopal, T., Serio, C., Melchionna, M., Esposito, A., Di
#'  Febbraro, M., & Raia, P. (2025). EutherianCoP. An integrated biotic and
#'  climate database for conservation paleobiology based on eutherian mammals.
#'  \emph{Scientific Data}, 12: 6. doi:10.1038/s41597-024-04181-4.
#'@examples
#' \dontrun{
#' library(RRgeo)
#'
#' setwd("YOUR_DIRECTORY")
#' getwd()->main.dir
#'
#' species_vec<-c( "Canis latrans","Canis lupus","Gulo gulo","Lutra lutra",
#'                 "Martes americana","Meles meles","Mustela erminea",
#'                 "Mustela nivalis","Procyon lotor","Ursus arctos","Ursus ingressus",
#'                 "Ursus maritimus","Ursus spelaeus","Vulpes velox","Vulpes vulpes" )
#'
#' eucop_data_preparation(input.dir=main.dir, species_name=species_vec,
#'                        variables="bio", calibration=TRUE, combine.ages="mean",
#'                        bk_points=list(),output.dir=main.dir)
#'
#'}



eucop_data_preparation<-function (input.dir,
                                  species_name,
                                  variables = "all",
                                  which.vars = NULL,
                                  calibration = FALSE,
                                  add.modern.occs = FALSE,
                                  combine.ages = NULL,
                                  remove.duplicates = TRUE,
                                  bk_points = NULL,
                                  output.dir = ".")
{
  message("Please cite:
          Mondanaro, A., Girardi, G., Castiglione, S., Timmermann, A., Zeller, E.,
          Venugopal, T., Serio, C., Melchionna, M., Esposito, A., Di Febbraro, M., &
          Raia, P. (2025). EutherianCoP. An integrated biotic and climate database
          for conservation paleobiology based on eutherian mammals.
          Scientific Data, 12: 6. doi:10.1038/s41597-024-04181-4.")
  misspacks <- sapply(c("rnaturalearth", "rnaturalearthdata",
                        "curl", "openxlsx"), requireNamespace, quietly = TRUE)
  if (any(!misspacks)) {
    stop("The following package/s are needed for this function to work, please install it/them:\n ",
         paste(names(misspacks)[which(!misspacks)], collapse = ", "),
         call. = FALSE)
  }
  if(!requireNamespace("rnaturalearthhires", quietly = TRUE)){
    stop("The package 'rnaturalearthhires' is needed for this function to work,\n please install it using the following command:
         install.packages('rnaturalearthhires', repos = 'http://packages.ropensci.org', type = 'source')")
  }
  if (calibration && (!requireNamespace("Bchron", quietly = TRUE)))
    stop("Package \"Bchron\" needed for calibration=TRUE. Please install it.",
         call. = FALSE)
  misspatstat <- sapply(c("spatstat.explore", "spatstat.geom"),
                        requireNamespace, quietly = TRUE)
  if (!is.null(bk_points) && any(!misspatstat))
    stop("Packages \"spatstat.geom\" and \"spatstat.explore\" needed to enable bk_points. Please install them.",
         call. = FALSE)
  biovar <- paste0("bio", c(1, 4:19))
  var_av <- c("LAI", "Megabiome", "NPP",
              paste0("PP", c("","JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG","SEP", "OCT", "NOV", "DEC")),
              paste0("TS",c("", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL","AUG", "SEP", "OCT", "NOV", "DEC")))
  if (variables == "bio") {
    if (!is.null(which.vars)) {
      if (all(!which.vars %in% biovar))
        stop(paste(which.vars, "is/are not available in the list of bioclimatic variables. Please see eucop_data_preparation documentation to ensure that the variable names match those available"))
      if (!all(which.vars %in% biovar)) {
        warning(paste(which.vars[!which.vars %in% biovar],
                      "is/are not available in the list of bioclimatic variables, it/they will be discarded"))
        which.vars <- which.vars[which.vars %in% biovar]
      }
    }else which.vars <- biovar
  }
  if (variables == "climveg") {
    if (!is.null(which.vars)) {
      if (all(!which.vars %in% var_av))
        stop(paste(which.vars, "is/are not available in the list of bioclimatic variables. Please see eucop_data_preparation documentation to ensure that the variable names match those available"))
      if (!all(which.vars %in% var_av)) {
        warning(paste(which.vars[!which.vars %in% var_av],
                      "is/are not available in the list of bioclimatic variables, it/they will be discarded"))
        which.vars <- which.vars[which.vars %in% var_av]
      }
    }else which.vars <- var_av
  }
  if (variables == "all") {
    if (!is.null(which.vars)) {
      if (all(!which.vars %in% c(var_av, biovar)))
        stop(paste(which.vars, "is/are not available in the list of bioclimatic variables. Please see eucop_data_preparation documentation to ensure that the variable names match those available"))
      if (!all(which.vars %in% c(var_av, biovar))) {
        warning(paste(which.vars[!which.vars %in% c(var_av,
                                                    biovar)], "is/are not available in the list of bioclimatic variables, it/they will be discarded"))
        which.vars <- which.vars[which.vars %in% c(var_av,
                                                   biovar)]
      }
    }else which.vars <- c(var_av, biovar)
  }

  curl::curl_download(url = "https://zenodo.org/records/14009105/files/1%20Raw%20occurrence%20data.xlsx?download=1",
                      destfile = paste0(input.dir, "/raw_data.xlsx"))
  occs <- openxlsx::read.xlsx(paste0(input.dir, "/raw_data.xlsx"),
                              sheet = 1)
  ages <- openxlsx::read.xlsx(paste0(input.dir, "/raw_data.xlsx"),
                              sheet = 2)
  occs_foss <- occs[occs$status == "fossil", ]
  spec <- split(occs_foss, occs_foss$spec)
  U <- unique(occs_foss$species)
  U <- sort(U)
  if (all(!species_name %in% U))
    stop("The selected species is/are not included in the fossil dataset")
  if (any(!species_name %in% U)) {
    missing <- species_name[!species_name %in% U]
    warning(paste("The selected species:", missing, "is/are not included in the fossil dataset"))
  }
  U <- U[U %in% species_name]
  occs_foss <- occs_foss[occs_foss$species %in% species_name,
  ]
  spec <- spec[names(spec) %in% species_name]
  if (variables == "climveg" | variables == "all") {
    if (file.exists(paste0(input.dir, "/climveg.zip")))
      print("climveg.zip file is already present in the input.dir folder") else
        curl::curl_download("https://zenodo.org/records/14936297/files/climveg.zip?download=1",
                            destfile = paste0(input.dir, "/climveg.zip"), quiet = FALSE)
    utils::unzip(paste0(input.dir, "/climveg.zip"), exdir = file.path(input.dir,"climveg"))
    nn <- paste0("X", seq(0, 130, 1), "kya")
    vv_bio <- list.files(paste0(input.dir, "/climveg"), full.names = TRUE,
                         recursive = TRUE)
    vv_bio <- vv_bio[grep(paste(nn, collapse = "|"), gsub(".tif",
                                                          "", basename(vv_bio)))]
    nn <- gsub(".tif", "", basename(vv_bio))
    vv_bio <- lapply(vv_bio, function(x) {
      xx <- rast(x)
      xx[[names(xx) %in% which.vars]]
    })
    names(vv_bio) <- nn
    vars <- vv_bio[mixedorder(names(vv_bio))]
  }
  if (variables == "bio" | variables == "all") {
    if (file.exists(paste0(input.dir, "/bio.zip")))
      print("bio.zip file is already present in the input.dir folder") else
        curl::curl_download("https://zenodo.org/records/14936297/files/bio.zip?download=1",
                            destfile = paste0(input.dir, "/bio.zip"), quiet = FALSE)
    utils::unzip(paste0(input.dir, "/bio.zip"), exdir = file.path(input.dir,"bio"))
    nn <- paste0("X", seq(0, 130, 1), "kya")
    vv_bio <- list.files(paste0(input.dir, "/bio"), full.names = TRUE,
                         recursive = TRUE)
    vv_bio <- vv_bio[grep(paste(nn, collapse = "|"), gsub(".tif",
                                                          "", basename(vv_bio)))]
    nn <- gsub(".tif", "", basename(vv_bio))
    vv_bio <- lapply(vv_bio, function(x) {
      xx <- rast(x)
      xx[[names(xx) %in% which.vars]]
    })
    names(vv_bio) <- nn
    vars1 <- vv_bio[mixedorder(names(vv_bio))]
    if (variables == "all") {
      vars <- lapply(1:length(vars1), function(x) {
        c(vars[[x]], vars1[[x]])
      })
      names(vars) <- names(vars1)
    }else vars <- vars1
  }

  dd <- occs_foss[match(ages$locality, occs_foss$locality),
  ]
  dd <- dd[-which(is.na(dd$locality)), ]
  datt <- ages[ages$locality %in% dd$locality, ]
  dat <- cbind(datt, longitude = dd$longitude, latitude = dd$latitude)
  dat$age <- as.numeric(dat$age)
  dat$uncertain <- as.numeric(dat$uncertain)
  dat$longitude <- round(as.numeric(as.character(dat$longitude)),
                         3)
  dat$latitude <- round(as.numeric(as.character(dat$latitude)),
                        3)
  dat_tot <- dat
  if (calibration) {
    data_tot <- dat_tot[grepl("AMS|14C", dat_tot$method),
    ]
    data_add <- dat_tot[!grepl("AMS|14C", dat_tot$method),
    ]
    data_add$curve <- "not calibrated"
    data_add$cal.age <- NA
    data_add$CI_2.5 <- NA
    data_add$CI_97.5 <- NA
    dating <- cal14C(data_tot, age = "age", uncertainty = "uncertain",
                         latitude = "latitude",
                         clust = NULL, save = FALSE)
    dating <- rbind(dating, data_add)
    colnames(dating)[colnames(dating) %in% c("cal.age", "CI_2.5",
                                             "CI_97.5")] <- c("mean.age", "min.age", "max.age")
    dating$mean.age[is.na(dating$mean.age)] <- dating$age[is.na(dating$mean.age)]
    dating$min.age[is.na(dating$min.age)] <- dating$age[is.na(dating$min.age)] -
      dating$uncertain[is.na(dating$min.age)]
    dating$max.age[is.na(dating$max.age)] <- dating$age[is.na(dating$max.age)] +
      dating$uncertain[is.na(dating$max.age)]
  }else {
    dating <- dat_tot
    dating$mean.age <- dating$age
    dating$min.age <- dating$uncertain
    dating$max.age <- dating$uncertain
  }
  datas <- lapply(U, function(x) dating[dating$locality %in%
                                          occs_foss[occs_foss$species == x, ][, "locality"], ])
  datas <- lapply(datas, function(x) x[order(x$locality), ])
  datas2 <- mapply(function(x, y) {
    x <- lapply(x$locality, function(z) {
      a <- do.call(rbind, replicate(nrow(y[y$locality ==
                                             z, ]), x[x$locality == z, ], simplify = FALSE))
      cbind(a, y[y$locality == z, ][, c("mean.age", "min.age",
                                        "max.age")])
    })
    x <- do.call(rbind, x)
  }, x = spec, y = datas, SIMPLIFY = FALSE)
  all_species <- do.call(rbind, datas2)
  all_sp_orig <- all_species[, c("species", "status", "longitude",
                                 "latitude", "locality", "country", "continent", "mean.age",
                                 "min.age", "max.age")]
  if (add.modern.occs) {
    occs_curr <- occs[occs$status == "current", ]
    if (species_name %in% unique(occs_curr$species)) {
      occs_curr <- occs_curr[which(occs_curr$species ==
                                     species_name), ]
      occs_curr$mean.age <- rep(0, nrow(occs_curr))
      occs_curr$min.age <- rep(0, nrow(occs_curr))
      occs_curr$max.age <- rep(0, nrow(occs_curr))
      occs_curr <- occs_curr[, c("species", "status", "longitude",
                                 "latitude", "locality", "country", "continent",
                                 "mean.age", "min.age", "max.age")]
      all_species <- rbind(all_sp_orig, occs_curr)
    }else all_species <- all_sp_orig
  }else all_species <- all_sp_orig
  if (!is.null(combine.ages)) {
    all_species <- split(all_species, all_species$species)
    all_species <- lapply(all_species, function(x) split(x,
                                                         x$locality))
    all_species <- lapply(all_species, function(x) do.call(rbind,
                                                           lapply(x, function(y) {
                                                             output <- y[1, ]
                                                             if (combine.ages == "mean")
                                                               output$mean.age <- round(mean(y$mean.age),
                                                                                        0)
                                                             if (combine.ages == "median")
                                                               output$mean.age <- round(median(y$mean.age),
                                                                                        0)
                                                             output
                                                           })))
    all_species <- do.call(rbind, all_species)
  }
  if (!is.null(bk_points)) {
    buff <- ifelse(is.null(bk_points$buff), 0.1, bk_points$buff)
    bk_strategy <- ifelse(is.null(bk_points$bk_strategy),
                          "background", bk_points$bk_strategy)
    bk_n <- ifelse(is.null(bk_points$bk_n), 10000, bk_points$bk_n)
    all_species <- all_species[, grep("max.age|min.age",
                                      colnames(all_species), invert = TRUE)]
    rownames(all_species) <- NULL
    sp <- split(all_species, all_species$species)
    lapply(1:length(sp), function(jj) {
      cat(paste("\n", "Formatting", unique(sp[[jj]]$species),
                "\n"))
      jjj2 <- st_as_sf(sp[[jj]], coords = c("longitude", "latitude"),
                       crs = 4326)
      jjj2 <- st_transform(jjj2,st_crs("ESRI:54009"))
      pol <- st_convex_hull(st_union(jjj2))
      buf <- max(st_distance(jjj2)) * buff
      pol <- st_buffer(pol, dist = as.numeric(buf))
      jjj2$time <- paste0("X", round(st_drop_geometry(jjj2)[, grepl("mean.age",
                                                                    colnames(jjj2))]/1000, 0), "kya")

      if (bk_strategy == "pseudoabsence") {
        mm <- vars[which(names(vars) %in% unique(jjj2$time))]
        mm <- sum(rast(lapply(mm, "[[", 1)), na.rm = TRUE)
        project(mm,st_crs(pol)$proj4string,res=50000)->mm
        mm <- mask(crop(mm, vect(pol)), vect(pol))
        dens <- density_background(jjj2, MASK = mm,rm.pres=TRUE)
      }
      ll <- split(jjj2, jjj2$time)
      rescale <- function(x, somma) {
        x * somma/sum(x)
      }
      WW0 <- sapply(ll, nrow)
      zeros <- ceiling(rescale(WW0, bk_n))
      if (sum(zeros)>bk_n){
        sum(zeros)-bk_n->rm
        zeros/sum(zeros)->pb
        sample(names(zeros), size = rm, replace = TRUE, prob = pb)->ppp
        for (group in unique(ppp)) {
          if (zeros[group] > 0) {
            zeros[group] <- zeros[group] - sum(ppp==group)
          }
        }
      }
      all <- lapply(ll, function(bb) {
        vv <- vars[[which(names(vars) %in% unique(bb$time))]]
        project(vv,st_crs(pol)$proj4string,res=50000)->vv
        vv <- mask(crop(vv, vect(pol)), vect(pol))
        xx <- fix.coastal.points(data = bb,
                                 r = vv[[1]], ncell = 2,
                                 occ.desaggregation = remove.duplicates)
        all <- data.frame(extract(vv, xx, ID = FALSE))
        all <- cbind(OBS = 1,xx, all)
        all <- all[, grep("time", colnames(all), invert = TRUE)]
        all <- all[!apply(is.na(all[, colnames(all) %in%
                                      names(vars[[1]]),drop=FALSE]), 1, all), ]
        all$mean.age <- round(all$mean.age/1000) * 1000

        if (bk_strategy == "pseudoabsence") {
          mm <- vv[[1]]
          mm[!is.na(mm), ] <- 1
          dens.ras <- dens * mm
          n.PA <- zeros[[which(names(zeros) %in% unique(bb$time))]]
          dens.ras.rtp <- cbind(data.frame(crds(dens.ras)),
                                as.points(dens.ras))
          dens.ras.rtp <- cbind(dens.ras.rtp, extract(dens.ras,
                                                      dens.ras.rtp[, 1:2], cells = TRUE)$cell)
          dens.ras.rtp[dens.ras.rtp[, 3] < 0.001, 3] <- 0.001
          PO <- xyFromCell(dens.ras, dens.ras.rtp[sample(1:nrow(dens.ras.rtp),
                                                         ifelse(nrow(dens.ras.rtp) < n.PA, nrow(dens.ras.rtp),
                                                                n.PA), prob = dens.ras.rtp[, 3]), 4])
          bk <- extract(vv, data.frame(PO), xy = TRUE,
                        ID = FALSE)
          st_as_sf(bk,coords=c("x","y"),crs=crs(mm))->vv1
        }

        if (bk_strategy == "background") {
          r1 <- extract(vv, all,
                        ID = FALSE, cells = TRUE, xy = TRUE)

          vv1 <- as.data.frame(vv, na.rm = TRUE, xy = TRUE,cells=TRUE)
          vv1[,match(colnames(r1),colnames(vv1))]->r2

          ad <- adjacent(vv[[1]], r1$cell, pairs = TRUE,
                         directions = 16)
          ad <- setdiff(unique(c(ad[, 2])), unique(c(ad[,
                                                        1])))
          ad <- r2[r2$cell %in% ad, ]
          ad <- ad[complete.cases(ad), , drop = FALSE]
          rr <- rbind(r1, ad)
          bk <- r2[!r2$cell %in% rr$cell, ]
          if (nrow(bk)>bk_n){
            bk <- bk[sample(nrow(bk), bk_n - nrow(rr)),
            ]
          }
          bk <- bk[, grep("cell", colnames(bk), invert = TRUE)]
          rr <- rr[, grep("cell", colnames(rr), invert = TRUE)]
          rbind(bk,rr)->bk
          st_as_sf(bk,coords=c("x","y"),crs=crs(vv))->vv1
        }
        vv1$OBS <- 0
        vv1$mean.age <- unique(all$mean.age)
        vv1$locality <- "bk_locs"
        vv1$continent <- NA
        vv1$country <- NA
        vv1$species <- unique(all$species)
        vv1$status <- names(which.max(table(all$status)))
        vv1 <- vv1[, match(colnames(all), colnames(vv1))]
        all <- rbind(all, vv1)
      })
      all_species <- do.call(rbind, all)
      rownames(all_species) <- NULL
      all_species[,colnames(all_species) %in% names(vars[[1]])] <-  apply(st_drop_geometry(all_species)[,
                                                                                                        colnames(all_species) %in% names(vars[[1]]),drop=FALSE],
                                                                          2, function(x) round(x, 3))


      colnames(all_species)[grepl("mean.age", colnames(all_species))] <- "age"
      suppressMessages(sf_use_s2(FALSE))
      ww <- rnaturalearth::ne_countries(returnclass = "sf", scale = 10)
      ww <- ww[!ww$continent == "Antarctica", ]
      st_transform(ww,st_crs("ESRI:54009"))->ww
      suppressWarnings(suppressMessages(ww <- st_crop(ww,
                                                      pol)))
      colnames(all_species)[colnames(all_species) == "continent"] <- "continent_old"
      all2 <- st_join(all_species, ww)[, c("geounit","continent")]
      all_species$country[!is.na(all2$geounit)] <- all2$geounit[!is.na(all2$geounit)]
      all_species$continent_old[!is.na(all2$continent)] <- all2$continent[!is.na(all2$continent)]
      colnames(all_species)[colnames(all_species) == "continent_old"] <- "continent"

      if (calibration && !is.null(combine.ages)) {
        if (file.exists(paste0(output.dir, "/", unique(all_species$species),
                               "_cal_combined.gpkg")))
          unlink(paste0(output.dir, "/", unique(all_species$species),
                        "_cal_combined.gpkg"))
        st_write(all_species, paste0(output.dir, "/",
                                     unique(all_species$species), "_cal_combined.gpkg"))
      }
      if (calibration && is.null(combine.ages)) {
        if (file.exists(paste0(output.dir, "/", unique(all_species$species),
                               "_cal_multi.gpkg"))) {
          unlink(paste0(output.dir, "/", unique(all_species$species),
                        "_cal_multi.gpkg"))
        }
        st_write(all_species, paste0(output.dir, "/",
                                     unique(all_species$species), "_cal_multi.gpkg"))
      }
      if (!calibration && !is.null(combine.ages)) {
        if (file.exists(paste0(output.dir, "/", unique(all_species$species),
                               "_uncal_combined.gpkg"))) {
          unlink(paste0(output.dir, "/", unique(all_species$species),
                        "_uncal_combined.gpkg"))
        }
        st_write(all_species, paste0(output.dir, "/",
                                     unique(all_species$species), "_uncal_combined.gpkg"))
      }
      if (!calibration && is.null(combine.ages)) {
        if (file.exists(paste0(output.dir, "/", unique(all_species$species),
                               "_uncal_multi.gpkg"))) {
          unlink(paste0(output.dir, "/", unique(all_species$species),
                        "_uncal_multi.gpkg"))
        }
        st_write(all_species, paste0(output.dir, "/",
                                     unique(all_species$species), "_uncal_multi.gpkg"))
      }
    })
  }else {
    rownames(all_species) <- NULL
    sp <- split(all_species, all_species$species)
    all_species <- lapply(1:length(sp), function(jj) {
      cat(paste("\n", "Formatting", unique(sp[[jj]]$species),
                "\n"))
      jjj2 <- st_as_sf(sp[[jj]], coords = c("longitude", "latitude"),
                       crs = 4326)
      jjj2 <- st_transform(jjj2,st_crs("ESRI:54009"))
      jjj2$time <- paste0("X", round(st_drop_geometry(jjj2)[, grepl("mean.age",
                                                                    colnames(jjj2))]/1000, 0), "kya")
      ll <- split(jjj2, jjj2$time)
      all <- lapply(ll, function(bb) {
        vv <- vars[[which(names(vars) %in% unique(bb$time))]]
        project(vv,st_crs(jjj2)$proj4string,res=50000)->vv
        xx <- fix.coastal.points(data = bb, r = vv, ncell = 2, occ.desaggregation = remove.duplicates)
        all <- extract(vv, xx, ID = FALSE)
        all <- cbind(OBS = 1, xx, all)
        all <- all[, grep("time", colnames(all), invert = TRUE)]
        all[!apply(is.na(all[, colnames(all) %in% names(vars[[1]]),drop=FALSE]),
                   1, all), ]
      })
      all_species <- do.call(rbind, all)
      rownames(all_species) <- NULL

      all_species[,colnames(all_species) %in% names(vars[[1]])] <-  apply(st_drop_geometry(all_species)[,
                                                                                                        colnames(all_species) %in% names(vars[[1]]),drop=FALSE],
                                                                          2, function(x) round(x, 3))

      colnames(all_species)[grepl("mean.age", colnames(all_species))] <- "age"
      suppressMessages(sf_use_s2(FALSE))
      ww <- rnaturalearth::ne_countries(returnclass = "sf", scale = 10)
      ww <- ww[!ww$continent == "Antarctica", ]
      st_transform(ww,st_crs("ESRI:54009"))->ww
      colnames(all_species)[colnames(all_species) == "continent"] <- "continent_old"
      all2 <- st_join(all_species, ww)[, c("geounit","continent")]
      all_species$country[!is.na(all2$geounit)] <- all2$geounit[!is.na(all2$geounit)]
      all_species$continent_old[!is.na(all2$continent)] <- all2$continent[!is.na(all2$continent)]
      colnames(all_species)[colnames(all_species) == "continent_old"] <- "continent"

      if (!calibration && is.null(combine.ages)) {
        all_species <- all_species[, grep("min.age",
                                          colnames(all_species), invert = TRUE)]
        colnames(all_species)[grepl("max.age", colnames(all_species))] <- "uncertain"
      }
      if (calibration && !is.null(combine.ages)) {
        all_species <- all_species[, grep("min.age|max.age",
                                          colnames(all_species), invert = TRUE)]
      }
      if (!calibration && is.null(combine.ages)) {
        all_species <- all_species[, grep("min.age|max.age",
                                          colnames(all_species), invert = TRUE)]
      }
      if (calibration && !is.null(combine.ages)) {
        if (file.exists(paste0(output.dir, "/", unique(all_species$species),
                               "_cal_combined.gpkg")))
          unlink(paste0(output.dir, "/", unique(all_species$species),
                        "_cal_combined.gpkg"))
        st_write(all_species, paste0(output.dir, "/",
                                     unique(all_species$species), "_cal_combined.gpkg"))
      }
      if (calibration && is.null(combine.ages)) {
        if (file.exists(paste0(output.dir, "/", unique(all_species$species),
                               "_cal_multi.gpkg"))) {
          unlink(paste0(output.dir, "/", unique(all_species$species),
                        "_cal_multi.gpkg"))
        }
        st_write(all_species, paste0(output.dir, "/",
                                     unique(all_species$species), "_cal_multi.gpkg"))
      }
      if (!calibration && !is.null(combine.ages)) {
        if (file.exists(paste0(output.dir, "/", unique(all_species$species),
                               "_uncal_combined.gpkg"))) {
          unlink(paste0(output.dir, "/", unique(all_species$species),
                        "_uncal_combined.gpkg"))
        }
        st_write(all_species, paste0(output.dir, "/",
                                     unique(all_species$species), "_uncal_combined.gpkg"))
      }
      if (!calibration && is.null(combine.ages)) {
        if (file.exists(paste0(output.dir, "/", unique(all_species$species),
                               "_uncal_multi.gpkg"))) {
          unlink(paste0(output.dir, "/", unique(all_species$species),
                        "_uncal_multi.gpkg"))
        }
        st_write(all_species, paste0(output.dir, "/",
                                     unique(all_species$species), "_uncal_multi.gpkg"))
      }
    })
  }
}

