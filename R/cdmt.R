#' Change Detection by Multispectral Trends
#'
#' This is the main function of the \pkg{cdmt} package.
#'
#' \code{cdmt} analyses inter-annual Landsat time series to detect changes in spectral trends at the pixel level.
#' Time series can be either univariate or multivariate, i.e. include multiple spectral bands/indices.
#' Impulsive noise, i.e. outliers in the time series, are removed through an iterative procedure.
#' One-year gaps in the time series are filled using either linear interpolation or extrapolation.
#' Input data are \code{SpatRaster} objects created by the \pkg{terra} package, which contain reflectance and/or spectral indices.
#' Changes in the intercept, slope or both of linear trends are detected using the High-dimensional Trend Segmentation (HiTS) procedure proposed by \insertCite{maeng2019adaptive;textual}{cdmt}.
#' The HiTS procedure aims at detecting changepoints in a piecewise linear signal where their number and location are unknown.
#'
#' @param sr_data character. Yearly surface reflectance data. Either the name of a \code{SpatRaster} object in the global environment or the full path where rasters (in \emph{.tif} format) are stored. Each layer name should include the year preceded by an underscore. If \code{NULL} only spectral indices are used.
#' @param si_data character. Yearly spectral indices. Either the name of a \code{SpatRaster} object in the global environment or the full path where rasters (in \emph{.tif} format) are stored. Each layer name should include the year preceded by an underscore. If \code{NULL} only reflectance bands are used.
#' @param out_path character. The path where output rasters (in \emph{.tif} format) are saved. Folders are created recursively if they do not exist. If \code{NULL} the output is a \code{SpatRaster} object created by the \pkg{terra} package.
#' @param sr numeric vector. Indices (positive integers) of the reflectance bands in each multiband raster to be included in the time series. If \code{NULL} all layers are used. Ignored if \code{sr_data} is a SpatRaster object.
#' @param si numeric vector. Indices (positive integers) of the spectral indices in each multiband raster to be included in the time series. If \code{NULL} all layers are used.
#' @param years numeric vector. Time interval (years) to be analysed.
#' @param cng_dir numeric vector. Direction of change caused by a disturbance in each spectral band/index. Valid values are either 1 or -1.
#' @param th_const numeric. Constant controlling the change threshold employed by the HiTS procedure during thresholding. Typical values are comprised in the interval \eqn{[1, 1.4]}.
#' @param noise_rm logical. If \code{TRUE} the impulsive noise filter is employed.
#' @param cores positive integer. Number of CPU cores used during analysis.
#'
#' @return A \code{SpatRaster} containing the following layers.
#'   \item{EST}{Estimated values (one layer for each year and band).}
#'   \item{CPT}{Detected changepoints (one layer for each year and band).}
#'   \item{SLO}{Slope of the linear segments (one layer for each year and band).}
#'   \item{MAG}{Magnitude in absolute terms (one layer for each year and band).}
#'   \item{MAG_REL}{Magnitude in relative terms (one layer for each year and band).}
#'   \item{LEN}{Length of segments (one layer per year).}
#'   \item{CPT_ID}{Type of change (one layer per year). One of the following values: 101 (abrupt disturbance); 102 (abrupt greening); 201 (gradual disturbance); 202 (gradual greening); 9 (other change).}
#'   \item{TPA}{Impulsive noise (one layer per year).}
#'   \item{D_DUR}{Duration of disturbance (one layer per year).}
#'   \item{G_DUR}{Duration of greening (one layer per year).}
#'   \item{D_MAX}{Maximum disturbance change magnitude (in relative tems) throughout the time series (one layer per band).}
#'   \item{D_FST}{Change magnitude (in relative tems) associated with the first disturbance detected within the time series (one layer per band).}
#'   \item{G_MAX}{Maximum greening change magnitude (in relative tems) throughout the time series (one layer per band).}
#'   \item{WGTS}{Weights assigned to each band (one layer per band).}
#'   \item{D_MAX_MD}{Median among bands using values of \code{D_MAX} (single layer).}
#'   \item{D_FST_MD}{Median among bands using values of \code{D_FST} (single layer).}
#'   \item{G_MAX_MD}{Median among bands using values of \code{G_MAX} (single layer).}
#'   \item{D_MAX_YR}{Year corresponding to \code{D_MAX_MD} (single layer).}
#'   \item{D_FST_YR}{Year corresponding to \code{D_FST_MD} (single layer).}
#'   \item{G_MAX_YR}{Year corresponding to \code{G_MAX_MD} (single layer).}
#'   \item{D_MAX_DR}{Duration of the disturbance corresponding to \code{D_MAX_MD} (single layer).}
#'   \item{D_FST_DR}{Duration of the disturbance corresponding to \code{D_FST_MD} (single layer).}
#'   \item{G_MAX_DR}{Duration of the greening corresponding to \code{G_MAX_MD} (single layer).}
#'   \item{D_MAX_ID}{Type of change (\code{CPT_ID}) of the disturbance corresponding to \code{D_MAX_MD} (single layer).}
#'   \item{D_FST_ID}{Type of change (\code{CPT_ID}) of the disturbance corresponding to \code{D_FST_MD} (single layer).}
#'   \item{G_MAX_ID}{Type of change (\code{CPT_ID}) of the greening corresponding to \code{G_MAX_MD} (single layer).}
#'   \item{N_GAP}{Number of gaps in the time series, if any (single layer).}
#'   \item{N_TPA}{Number of years containing impulsive noise, if any (single layer).}
#'   \item{N_ITER}{Number of iterations performed by the impulsive noise filter (single layer).}
#' If a valid \code{out_path} is provided, output rasters in \emph{*.tif} format are directly written in the output folder.
#'
#' @author Donato Morresi, \email{donato.morresi@@gmail.com}
#'
#' @seealso \code{\link{cdmt_int}}
#'
#' @references
#' \insertRef{maeng2019adaptive}{cdmt}
#'
#' @examples
#' # Load raster data
#' data(lnd_sr)
#' data(lnd_si)
#' lnd_sr <- terra::rast(lnd_sr)
#' lnd_si <- terra::rast(lnd_si)
#'
#' # Process data
#' rsout <- cdmt(sr_data = "lnd_sr",
#'               si_data = "lnd_si",
#'               years = 1985:2020,
#'               cng_dir = c(-1, -1, 1, 1),
#'               cores = 2
#'               )
#'
#' # Plot the median spectral change magnitude across bands (D_MAX_MD)
#' terra::plot(rsout[["D_MAX_MD"]])
#'
#' @export
#' @import matrixStats
#' @import parallel
#' @importFrom accelerometry rle2
#' @importFrom lubridate seconds_to_period
#' @importFrom Rdpack reprompt
#' @importFrom stats .lm.fit
#' @importFrom stats median
#' @importFrom terra app
#' @importFrom terra nlyr
#' @importFrom terra rast
#' @importFrom terra tmpFiles
#' @importFrom terra writeRaster


cdmt <- function(sr_data, si_data, out_path = NULL, sr = NULL, si = NULL, years, cng_dir, th_const = 1.2, noise_rm = TRUE, cores = 1) {

  message("\nLoading raster data ...")

  if (is.null(sr_data) && is.null(si_data)) {
    stop("missing input data")
  }

  if (!is.null(sr_data)) {

    if (exists(sr_data, envir = .GlobalEnv)) {  #&& inherits(sr_data, "SpatRaster")
      sr_rs <- get(sr_data, envir = .GlobalEnv)

      if (!inherits(sr_rs, "SpatRaster")) {
        stop("unrecognised reflectance data")
      }
    }
    else if (inherits(sr_data, "connection")) {
      pat <- paste0("(", paste0("_", years, collapse = "|"), ")+.+tif$")
      f <- list.files(sr_data, pat, full.names = TRUE)
      sr_rs <- lapply(f, function(x) rast(x, lyrs = sr))
      sr_rs <- do.call(c, sr_rs)
    }
    else {
      stop("unrecognised reflectance data")
    }
  }

  if (!is.null(si_data)) {

    if (exists(si_data, envir = .GlobalEnv)) {
      si_rs <- get(si_data, envir = .GlobalEnv)

      if (!inherits(si_rs, "SpatRaster")) {
        stop("unrecognised spectral indices data")
      }
    }
    else if (inherits(si_data, "connection")) {
      pat <- paste0("(", paste0("_", years, collapse = "|"), ")+.+tif$")
      f <- list.files(si_data, pat, full.names = TRUE)
      si_rs <- lapply(f, function(x) rast(x, lyrs = si))
      si_rs <- do.call(c, si_rs)
    }
    else {
      stop("unrecognised spectral indices data")
    }
  }

  if (!is.null(sr_data) && !is.null(si_data)) {

    rs <- lapply(years, function(y) {
      sr_y <- sr_rs[[grep(paste0("_", y), names(sr_rs))]]
      si_y <- si_rs[[grep(paste0("_", y), names(si_rs))]]
      rs_y <- c(sr_y, si_y)
    })

    rs <- do.call(c, rs)
  }
  else if (!is.null(sr_data)) {
    rs <- sr_rs
  }
  else {
    rs <- si_rs
  }

  # Output layer names
  matn <- c("EST", "CPT", "SLO", "MAG", "MAG_REL")               # matrices
  ve1n <- c("LEN", "CPT_ID", "TPA", "D_DUR", "G_DUR")            # long vectors (years)
  ve2n <- c("D_MAX", "D_FST", "G_MAX", "WGTS")                   # short vectors (bands)
  valn <- c("D_MAX_MD", "D_FST_MD", "G_MAX_MD", "D_MAX_YR",      # single values
            "D_FST_YR", "G_MAX_YR", "D_MAX_DR", "D_FST_DR",
            "G_MAX_DR", "D_MAX_ID", "D_FST_ID", "G_MAX_ID",
            "N_GAP", "N_TPA", "N_ITER")

  ny <- length(years)
  nb <- nlyr(rs) / ny

  if (nb %% 1 != 0) {
    stop("unequal number of raster layers")
  }

  if (length(cng_dir) != nb) {
    stop("length of cng_dir does not match the number of bands")
  }

  # Compute number of output layers
  nout <- (length(matn) * nb * ny + length(ve1n) * ny + length(ve2n) * nb + length(valn))

  # Create empty vector for invalid pixels
  ev <- rep(NA, nout)

  message("\nProcessing Landsat time series ...")

  if (cores > 1L) {

    cl <- makeCluster(cores)
    on.exit(stopCluster(cl), add=TRUE)

    clusterEvalQ(cl, {
      library(cdmt)
    })

    tic <- proc.time()
    rsout <- app(rs,
                 fun = cdmt_int,
                 nb = nb,
                 ny = ny,
                 yrs = years,
                 ev = ev,
                 dir = cng_dir,
                 th_const = th_const,
                 noise_rm = noise_rm,
                 cores = cl)
    toc <- proc.time()
  }
  else {

    tic <- proc.time()
    rsout <-  app(rs,
                  fun = cdmt_int,
                  nb = nb,
                  ny = ny,
                  yrs = years,
                  ev = ev,
                  dir = cng_dir,
                  th_const = th_const,
                  noise_rm = noise_rm)
    toc <- proc.time()
  }

  elapsed <- seconds_to_period((toc - tic)[3])

  # Compute indices and set layer names
  ni <- seq_len(length(rep(seq_len(nb), ny)) * length(matn))
  mi <- matrix(ni, ncol = nb, byrow = TRUE)
  mn <- paste0(rep(matn, each=ny), "_", years)

  v1i <- (ni[length(ni)] + 1):(ni[length(ni)] + ny * length(ve1n))
  v1im <- matrix(v1i, ncol = ny, byrow = TRUE)
  v1n <- c(paste0(ve1n, "_", years[1], "_", years[ny]))

  v2i <- (v1i[length(v1i)] + 1):(v1i[length(v1i)] + nb * length(ve2n))
  v2im <- matrix(v2i, ncol = nb, byrow = TRUE)
  v2n <- c(paste0(ve2n, "_", years[1], "_", years[ny]))

  si <- (v2i[length(v2i)] + 1):(v2i[length(v2i)] + length(valn))
  sim <- matrix(si, ncol = 1)
  sn <- paste0(valn)

  out_ind <- list(mi, v1im, v2im, sim)
  out_nm <- list(mn, v1n, v2n, sn)

  names(rsout) <- c(paste0(rep(mn, each=nb), "_", seq(nb)),
                    paste0(rep(v1n, each=ny), "_", seq(ny)),
                    paste0(rep(v2n, each=nb), "_", seq(nb)),
                    sn)

  if (is.null(out_path)) {

    message("\nExecution took ", elapsed)
    tmpFiles(remove = TRUE)
    return(rsout)
  }
  else {

    if (!inherits(out_path, "connection")) {
      stop("invalid out_path")
    }

    if (!dir.exists(out_path)) {
      dir.create(out_path, recursive = TRUE)
    }

    message("\nWrite output rasters to disk ...")

    for (i in seq_along(out_ind)) {

      ind <- out_ind[[i]]
      nm <- out_nm[[i]]

      for (j in seq_along(nm)) {
        rs <- rsout[[ind[j,]]]
        writeRaster(rs, file.path(paste0(out_path, "/", nm[j], ".tif")), overwrite = TRUE,
                    wopt = list(filetype = "GTiff", datatype = "FLT4S", NAflag = -32768, gdal = c("COMPRESS=LZW")))
      }
    }
    message("\nExecution took ", elapsed)
    tmpFiles(remove = TRUE)
  }
}
