#' Change Detection by Multispectral Trends
#'
#' This is the main function of the package.
#'
#' \code{cdmt()} detects changes in linear trends of inter-annual Landsat time series at the pixel level.
#' Time series can be either univariate or multivariate, \emph{i.e.} include multiple spectral bands/indices.
#' Changes in the intercept, slope or both of linear trends are detected using the High-dimensional Trend Segmentation (HiTS) procedure proposed by Maeng (2019).
#' The HiTS procedure aims at detecting changepoints in a piecewise linear signal where the number and location of changepoints are unknown.
#' Input and output raster data is handled using the \pkg{terra} package.
#'
#' @param sr_data character. Yearly surface reflectance data. Either the name of a \code{SpatRaster} object in the global environment or the full path where rasters (in \emph{*.tif} format) are stored. Each layer name should include the year preceded by an underscore. If \code{NULL} only spectral indices are used.
#' @param si_data character. Yearly spectral indices. Either the name of a \code{SpatRaster} object in the global environment or the full path where rasters (in \emph{*.tif} format) are stored. Each layer name should include the year preceded by an underscore. If \code{NULL} only reflectance bands are used.
#' @param out_path character. The path where output rasters (in \emph{*.tif} format) are saved. Folders are created recursively if they do not exist. If \code{NULL} the output is a \code{SpatRaster} object created by the \pkg{terra} package.
#' @param sr numeric vector. Indices (positive integers) of the reflectance bands in each multiband raster to be included in the time series. If \code{NULL} all layers are used. Ignored if \code{sr_data} is a SpatRaster object.
#' @param si numeric vector. Indices (positive integers) of the spectral indices in each multiband raster to be included in the time series. If \code{NULL} all layers are used.
#' @param years numeric vector. Time interval (years) to be analysed.
#' @param cng_dir numeric vector. Direction of change caused by a disturbance in each spectral band/index. Valid values are either 1 or -1.
#' @param th_const numeric. Constant controlling the change threshold employed by the HiTS procedure during thresholding. Typical values are comprised in the interval \eqn{[1, 1.4]}.
#' @param noise_rm logical. If \code{TRUE} the impulsive noise filter is active.
#' @param cores positive integer. Number of CPU cores used during analysis.
#'
#' @return A \code{SpatRaster} object containing the following layers.
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
#' If a valid \code{out_path} is provided, a series of rasters in \emph{*.tif} format containing either one or multiple layers are directly written in the output folder.
#'
#' @author Donato Morresi, \email{donato.morresi@@gmail.com}
#'
#' @export
#' @import matrixStats
#' @import parallel
#' @importFrom accelerometry rle2
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
    stop('\nExecution halted: missing raster data')
  }

  if (!is.null(sr_data)) {

    if (exists(sr_data, .GlobalEnv)) {
      sr_rs <- get(sr_data, envir = .GlobalEnv)
    }
    else {
      pat <- paste0("(", paste0("_", years, collapse = "|"), ")+.+tif$")
      f <- list.files(sr_data, pat, full.names = TRUE)
      sr_rs <- lapply(f, function(x) rast(x, lyrs = sr))
      sr_rs <- do.call(c, sr_rs)
    }
  }

  if (!is.null(si_data)) {

    if (exists(si_data, .GlobalEnv)) {
      si_rs <- get(si_data, envir = .GlobalEnv)
    }
    else {
      pat <- paste0("(", paste0("_", years, collapse = "|"), ")+.+tif$")
      f <- list.files(si_data, pat, full.names = TRUE)
      si_rs <- lapply(f, function(x) rast(x, lyrs = si))
      si_rs <- do.call(c, si_rs)
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
    stop("Execution halted: unequal number of raster layers")
  }

  if (length(cng_dir) != nb) {
    stop("Execution halted: the length of cng_dir does not match the number of bands")
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

  elapsed <- round((toc - tic)[3]/3600, 2)

  # Compute indices and set layer names
  ni <- seq_len(length(rep(seq_len(nb), ny)) * length(matn))
  mi <- matrix(ni, ncol = nb, byrow = TRUE)
  mn <- paste0("CDMT_", years, "_", rep(matn, each=ny))

  v1i <- (ni[length(ni)]+1):(ni[length(ni)] + ny * length(ve1n))
  v1im <- matrix(v1i, ncol = ny, byrow = TRUE)
  v1n <- c(paste0("CDMT_", ve1n, "_", years[1], "_", years[ny]))

  v2i <- (v1i[length(v1i)] + 1):(v1i[length(v1i)] + nb * length(ve2n))
  v2im <- matrix(v2i, ncol = nb, byrow = TRUE)
  v2n <- c(paste0("CDMT_", ve2n, "_", years[1], "_", years[ny]))

  si <- (v2i[length(v2i)] + 1):(v2i[length(v2i)] + length(valn))
  sn <- paste0("CDMT_", valn)

  names(rsout) <- c(paste0(rep(mn, each=nb), "_", seq(nb)),
                    paste0(rep(v1n, each=ny), "_", seq(ny)),
                    paste0(rep(v2n, each=nb), "_", seq(nb)),
                    sn)

  if (!is.null(out_path)) {

    if (!dir.exists(out_path)) {
      dir.create(out_path, recursive = TRUE)
    }

    message("\nWrite output rasters to disk ...")

    # Write data stored as matrices
    lapply(seq_along(mn), function(i) {
      rs <- rsout[[mi[i,]]]
      writeRaster(rs, file.path(paste0(out_path, "/", mn[i], ".tif")), overwrite = TRUE,
                  wopt = list(filetype = "GTiff", datatype = "FLT4S", NAflag = -32768, gdal = c("COMPRESS=LZW")))
    })

    # Write data stored as long vectors
    lapply(seq_along(v1n), function(i) {
      rs <- rsout[[v1im[i,]]]
      writeRaster(rs, file.path(paste0(out_path, "/", v1n[i], ".tif")), overwrite = TRUE,
                  wopt = list(filetype = "GTiff", datatype = "FLT4S", NAflag = -32768, gdal = c("COMPRESS=LZW")))
    })

    # Write data stored as short vectors
    lapply(seq_along(v2n), function(i) {
      rs <- rsout[[v2im[i,]]]
      writeRaster(rs, file.path(paste0(out_path, "/", v2n[i], ".tif")), overwrite = TRUE,
                  wopt = list(filetype = "GTiff", datatype = "FLT4S", NAflag = -32768, gdal = c("COMPRESS=LZW")))
    })

    # Write data stored as single values
    lapply(seq_along(sn), function(i) {
      rs <- rsout[[si[i]]]
      writeRaster(rs, file.path(paste0(out_path, "/", sn[i], ".tif")), overwrite = TRUE,
                  wopt = list(filetype = "GTiff", datatype = "FLT4S", NAflag = -32768, gdal = c("COMPRESS=LZW")))
    })
  }
  else {
    return(rsout)
  }

  message("\nExecution took ", elapsed, " hours")

  tmpFiles(remove = TRUE)
}
