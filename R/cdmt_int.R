#' Pixel-wise Change Detection by Multispectral Trends
#'
#' This is the internal function called by \code{\link{cdmt}}.
#'
#' \code{cdmt_int} is not directly called by the user.
#' It analyses inter-annual Landsat time series to detect changes in spectral trends at the pixel level.
#' Time series can be either univariate or multivariate, i.e. include multiple spectral bands/indices.
#' Input data can be a \code{numeric} \code{vector} extracted from a \code{SpatRaster} object or a \code{matrix}.
#' Impulsive noise, i.e. outliers in the time series, are removed through an iterative procedure.
#' One-year gaps in the time series are filled using either linear interpolation or extrapolation.
#' Changes in the intercept, slope or both of linear trends are detected using the High-dimensional Trend Segmentation (HiTS) procedure proposed by \insertCite{maeng2019adaptive;textual}{cdmt}.
#' The HiTS procedure aims at detecting changepoints in a piecewise linear signal where their number and location are unknown.
#'
#' @param x matrix. Each row contains time-ordered data relative to an individual reflectance band or spectral index. Each column contains data relative to a single year.
#' @param nb integer. The number of reflectance bands and/or spectral indices, i.e. rows of \code{x}.
#' @param ny integer. The number of years, i.e. columns of \code{x}.
#' @param yrs numeric vector. The sequence of years to be analysed.
#' @param ev logical vector. A vector containing \code{NA} values to be used in case of processing failure.
#' @param dir numeric vector. Direction of change caused by a disturbance in each spectral band/index. Valid values are either 1 or -1.
#' @param th_const numeric. Constant controlling the change threshold employed by the HiTS procedure during thresholding. Typical values are comprised in the interval \eqn{[1, 1.4]}.
#' @param noise_rm logical. If \code{TRUE} the impulsive noise filter is used.
#' @param as_list logical. If If \code{TRUE} the output is a \code{list}. Otherwise, the output is a \code{numeric} \code{vector}.
#'
#' @return If \code{as_list} is \code{TRUE}, a \code{list} containing the following elements.
#'   \item{est}{Estimated values (one layer for each year and band).}
#'   \item{cpt}{Detected changepoints (one layer for each year and band).}
#'   \item{slo}{Slope of the linear segments (one layer for each year and band).}
#'   \item{mag}{Magnitude in absolute terms (one layer for each year and band).}
#'   \item{mag_rel}{Magnitude in relative terms (one layer for each year and band).}
#'   \item{len}{Length of segments (one layer per year).}
#'   \item{cpt_id}{Type of change (one layer per year). One of the following values: 101 (abrupt disturbance); 102 (abrupt greening); 201 (gradual disturbance); 202 (gradual greening); 9 (other change).}
#'   \item{tpa}{Impulsive noise (one layer per year).}
#'   \item{d_dur}{Duration of disturbance (one layer per year).}
#'   \item{g_dur}{Duration of greening (one layer per year).}
#'   \item{d_max}{Maximum disturbance change magnitude (in relative tems) throughout the time series (one layer per band).}
#'   \item{d_fst}{Change magnitude (in relative tems) associated with the first disturbance detected within the time series (one layer per band).}
#'   \item{g_max}{Maximum greening change magnitude (in relative tems) throughout the time series (one layer per band).}
#'   \item{wgts}{Weights assigned to each band (one layer per band).}
#'   \item{d_max_md}{Median among bands using values of \code{D_MAX} (single layer).}
#'   \item{d_fst_md}{Median among bands using values of \code{D_FST} (single layer).}
#'   \item{g_max_md}{Median among bands using values of \code{G_MAX} (single layer).}
#'   \item{d_max_yr}{Year corresponding to \code{D_MAX_MD} (single layer).}
#'   \item{d_fst_yr}{Year corresponding to \code{D_FST_MD} (single layer).}
#'   \item{g_max_yr}{Year corresponding to \code{G_MAX_MD} (single layer).}
#'   \item{d_max_dr}{Duration of the disturbance corresponding to \code{D_MAX_MD} (single layer).}
#'   \item{d_fst_dr}{Duration of the disturbance corresponding to \code{D_FST_MD} (single layer).}
#'   \item{g_max_dr}{Duration of the greening corresponding to \code{G_MAX_MD} (single layer).}
#'   \item{d_max-id}{Type of change (\code{CPT_ID}) of the disturbance corresponding to \code{D_MAX_MD} (single layer).}
#'   \item{d_fst_id}{Type of change (\code{CPT_ID}) of the disturbance corresponding to \code{D_FST_MD} (single layer).}
#'   \item{g_max_id}{Type of change (\code{CPT_ID}) of the greening corresponding to \code{G_MAX_MD} (single layer).}
#'   \item{n_gap}{Number of gaps in the time series, if any (single layer).}
#'   \item{n_tpa}{Number of years containing impulsive noise, if any (single layer).}
#'   \item{n_iter}{Number of iterations performed by the impulsive noise filter (single layer).}
#' If \code{as_list} is \code{FALSE} the result is a \code{numeric} \code{vector}.
#'
#' @author Donato Morresi, \email{donato.morresi@@gmail.com}
#'
#' @references
#' \insertRef{maeng2019adaptive}{cdmt}
#'
#' @examples
#' # Load raster data
#' data(lnd_si)
#' lnd_si <- terra::rast(lnd_si)
#'
#' # Extract values
#' v <- terra::values(lnd_si)[2000,]
#'
#' # Process data
#' out <- cdmt_int(v,
#'                 nb = 3,
#'                 ny = 36,
#'                 yrs = 1985:2020,
#'                 dir = c(-1, 1, 1),
#'                 as_list = TRUE
#'                 )
#'
#' # Plot results
#' var <- c("MSI", "TCW", "TCA")
#' par(mfrow=c(length(var), 1))
#'
#' for(i in seq_along(var)) {
#'   plot(matrix(v, nrow = 3, ncol = 36)[i,], type = "l", col = 2,
#'   lwd = 2, ylab = "", xlab = "", xaxt = "n")
#'   lines(out$est[i,], type = "l", col = 4, lwd = 2, lty = 1)
#'   abline(v = which(out$cpt_id == 101), lty = 2, col = 1, lwd = 2)
#'   axis(1, at = seq_along(1985:2020), labels = 1985:2020)
#'   title(main = paste(var[i]), adj = 0)
#' }
#'
#' @export

cdmt_int <- function(x, nb, ny, yrs, ev = NA, dir, th_const = 1.2, noise_rm = TRUE, as_list = FALSE) {

  if (all(is.na(x))) {
    return(ev)
  }

  dim(x) <- c(nb, ny)

  cc <- !colAnyNAs(x)

  ncc <- sum(cc)

  if (ncc < ny) {

    gl <- rle2(cc, class='logical')
    gl <- gl[gl[,1]==0L, 2]

    if (any(gl > 1L)) {
      return(ev)
    }
    else {
      x <- x[, cc, drop=F]
    }
  }

  sd <- sd_est(x)

  if (any(sd == 0L)) {
    return(ev)
  }

  w <- wr2(x)

  tpa.iter <- 1L
  ntpa <- iter <- 0L
  tpav <- integer(dim(x)[2])

  if (noise_rm) {

    while (tpa.iter > 0L & iter < 6L) {

      iter <- iter + 1L
      tpa.iter <- 0L

      if (any(tpav != 0L)) {

        # estimate standard deviation without tpas
        sd <- sd_est(x[, tpav == 0L, drop=F])

        if (any(sd == 0L)) {
          return(ev)
        }

        # compute weights
        w <- wr2(x[, tpav == 0L, drop=F])
      }

      # perform hits
      y <- hd.bts.cpt(x, sd, th_const, w)
      tpa <- 0L

      # check for the presence of noise
      if (y$no.of.cpt > 0L) {

        cpt <- cpt_cnd(y, x, sd)

        # check PAs if present among change points
        if (length(cpt) > 0L) {

          tpa <- proc_cpt(x, sd, y, cpt, th_const, w, dir)
        }

        if (any(tpa > 0L)) {

          tpa.iter <- length(tpa)
          ntpa <- sum(ntpa, tpa.iter)
          tpav[tpa] <- 1L

          for (i in seq_along(tpa)) {

            tpai <- tpa[i]
            x[,tpai] <- colMeans2(rbind(x[,tpai - 1L], x[,tpai + 1L]))
          }
        }
      }
    }

  }
  else {
    y <- hd.bts.cpt(x, sd, th_const, w)
  }

  if (ncc < ny) {

    ec <- which(cc==F)

    # update cpts if there were data gaps
    if (y$no.of.cpt > 0L) {

      fc <- seq_len(ny)[cc]

      for (i in seq_along(y$cpt)) {

        y$cpt[i] <- fc[y$cpt[i]]
      }
    }

    # fill missing columns
    for (i in seq_along(ec)) {

      eci <- ec[i]

      # process gaps occurring at the beginning
      if (eci == 1L) {

        x <- cbind(rep(NA, nb), x)
        y$est <- cbind(rep(NA, nb), y$est)
        tpav <- c(0L, tpav)

        if (any(y$cpt != 0L) && any(y$cpt == eci + 1L || y$cpt == eci + 2L)) {
          # use the following value
          x[,eci] <- y$est[,eci] <- y$est[,eci + 1L]
        }
        else {
          # use extrapolated values
          x[,eci] <- y$est[,eci] <- y$est[,eci + 2L] + 2*(y$est[,eci + 1L] - y$est[,eci + 2L])
        }
      }
      # or at the end
      else if (eci == ny) {

        x <- cbind(x, rep(NA, nb))
        y$est <- cbind(y$est, rep(NA, nb))
        tpav <- c(tpav, 0L)

        if (any(y$cpt != 0L) && any(y$cpt == eci - 1L || y$cpt == eci - 2L)) {
          # use the preceding value
          x[,eci] <- y$est[,eci] <- y$est[,eci - 1L]
        }
        else {
          # use extrapolated values
          x[,eci] <- y$est[,eci] <- y$est[,eci - 2L] + 2*(y$est[,eci - 1L] - y$est[,eci - 2L])
        }
      }
      else {
        seq1 <- seq_len(eci - 1L)
        seq2 <- eci:dim(x)[2]
        x <- cbind(x[,seq1,drop=F], rep(NA, nb), x[,seq2,drop=F])
        y$est <- cbind(y$est[,seq1,drop=F], rep(NA, nb), y$est[,seq2,drop=F])
        tpav <- c(tpav[seq1], 0L, tpav[seq2])

        if (any(y$cpt != 0L) && any(y$cpt == eci + 1L)) {

          if (eci > 2L) {

            # fill gap using values extrapolated from the preceding segment
            x[,eci] <- y$est[,eci] <- y$est[,eci - 2L] + 2*(y$est[,eci - 1L] - y$est[,eci - 2L])
          }
          else {
            # use the preceding value
            x[,eci] <- y$est[,eci] <- y$est[,eci - 1L]
          }
        }
        else {

          # perform linear interpolation if there aren't PA or CPT near the gap
          x[,eci] <- y$est[,eci] <- colMeans2(rbind(y$est[,eci - 1L], y$est[,eci + 1L]))
        }
      }
    }
  }

  # compute the number of segments
  nseg <- y$no.of.cpt + 1L

  # compute the number of gaps
  ngap <- ny - ncc

  # store cptind per each year in a matrix
  cpts <- matrix(NA, nb, ny)
  #cpts[,y$cpt] <- y$cptind

  # assign the year to the gap occurrence
  tpav <- tpav * yrs

  #strtoi(paste(cpts, collapse=''), base=2)

  # compute length of segments
  if (nseg == 1L) {
    len <- ny
  }
  else {
    len <- integer(nseg)

    for (i in seq_len(nseg)) {
      if (i == 1L) {
        len[i] <- y$cpt[1] - 1L
      }
      else if (i == nseg) {
        len[i] <- ny - y$cpt[y$no.of.cpt] + 1L
      }
      else {
        len[i] <- y$cpt[i] - y$cpt[i - 1L]
      }
    }
  }

  # create vector with (repeated) length of segments
  len.seg <- rep(len, times=len)

  # compute the slope of each segment
  slo <- matrix(NA, nb, nseg)

  if (nseg == 1L) {
    slo[,1] <- y$est[,2] - y$est[,1]
  }
  else {
    for (i in seq_len(nseg)) {

      if (len[i] == 1L) {
        slo[,i] <- rep(0, nb)
      }
      else if (i == nseg) {
        slo[,i] = y$est[,ny] - y$est[,ny-1L]
      }
      else {
        slo[,i] <- y$est[,y$cpt[i]-1L] - y$est[,y$cpt[i]-2L]
      }
    }
  }

  # repeat columns based on the length of each segment
  slo.seg <- matrix(NA, nb, ny)

  for (i in seq_len(nb)) {
    slo.seg[i,] <- rep(slo[i,], len)
  }

  # compute the magnitude of change and identify the type of change
  mag <- mag.rel <- matrix(NA, nb, ny)
  cpt.id <- rep(NA, ny)
  d.mag.co <- d.fst.mg <- g.mag.co <- rep(NA, nb)
  d.mag.yr <- d.fst.yr <- g.mag.yr <- NA
  d.mag <- d.fst.md <- g.mag <- NA

  d.dur <- g.dur <- rep(NA, ny)
  d.mag.dr <- g.mag.dr <- d.fst.dr <- d.mag.id <- d.fst.id <- g.mag.id <- NA

  if (nseg > 1L) {

    for (i in seq_along(y$cpt)) {

      cpti <- y$cpt[i]

      # ignore changes occurring at the end of the time series
      if (cpti == ny) {
        next
      }

      mag1.act <- x[,cpti-1L] - x[,cpti]
      mag1.est <- y$est[,cpti-1L] - y$est[,cpti]

      mag2.act <- x[,cpti-1L] - x[,cpti+1L]
      mag2.est <- y$est[,cpti-1L] - y$est[,cpti+1L]

      cng.dir.act <- sign(mag1.act)
      cng.dir.est <- sign(mag1.est)

      ### assess if the change is abrupt or gradual
      if (len.seg[cpti] == 1L) {
        cng.type <- 1L  # abrupt change
      }
      else if (sum(abs(mag1.est) > sd, na.rm=T) > (nb/2) && sum(abs(mag2.est) > sd, na.rm=T) > (nb/2) && sum(sign(mag1.est) == sign(mag2.est)) > (nb/2)) {
        cng.type <- 1L  # abrupt change
      }
      else {
        cng.type <- 2L   # gradual change
      }

      if (cng.type == 1L) {   # abrupt change

        if (sum(cng.dir.act == dir) > nb/2 && sum(cng.dir.est == dir) > nb/2) {    # pattern of change corresponding to a disturbance

          cpt.id[cpti] <- 101L    # abrupt disturbance
          d.dur[cpti] <- 1L
          mag[,cpti] <- mag1.est
          mag.rel[,cpti] <- mag1.est/y$est[,cpti-1L] * 100

        }
        else if (sum(cng.dir.act == -dir) > nb/2 && sum(cng.dir.est == -dir) > nb/2) {

          cpt.id[cpti] <- 102L    # abrupt greening
          g.dur[cpti] <- len.seg[cpti]
          mag[,cpti] <- mag1.est
          mag.rel[,cpti] <- mag1.est/y$est[,cpti-1L] * 100

        }
        else {
          cpt.id[cpti] <- 9L
        }
      }
      else {   # gradual change

        slo.dir <- sign(slo.seg[,cpti])
        mag3.est <- y$est[,cpti] - y$est[,cpti-1L + len.seg[cpti]]

        if (sum(slo.dir == -dir, na.rm=T) > (nb/2)) {
          cpt.id[cpti] <- 201L    # gradual decline
          d.dur[cpti] <- len.seg[cpti]
          mag[,cpti] <- mag3.est
          mag.rel[,cpti] <- mag3.est/y$est[,cpti-1L] * 100
        }
        else if (sum(slo.dir == dir, na.rm=T) > (nb/2)) {
          cpt.id[cpti] <- 202L    # gradual greening
          g.dur[cpti] <- len.seg[cpti]
          mag[,cpti] <- mag3.est
          mag.rel[,cpti] <- mag3.est/y$est[,cpti-1L] * 100
        }
        else {
          cpt.id[cpti] <- 9L
        }
      }
    }

    # find the highest-magnitude disturbance and its duration
    mag.md <- colMedians(abs(mag.rel))

    if (any(cpt.id %in% c(101, 201))) {
      d.mag <- max(mag.md[cpt.id %in% c(101, 201)])
      d.mag.ci <- which(mag.md %in% d.mag)
      d.mag.co <- mag.rel[, d.mag.ci]
      d.mag.yr <- yrs[d.mag.ci]
      d.mag.dr <- d.dur[d.mag.ci]
      d.mag.id <- cpt.id[d.mag.ci]

      # find the first disturbance occurred
      ind <- which.max(cpt.id %in% c(101, 201))
      d.fst.yr <- yrs[ind]
      d.fst.mg <- mag.rel[, ind]
      d.fst.md <- median(d.fst.mg)
      d.fst.dr <- d.dur[ind]
      d.fst.id <- cpt.id[ind]
    }

    # find the highest-magnitude greening and its duration
    if (any(cpt.id %in% c(102, 202))) {
      g.mag <- max(mag.md[cpt.id %in% c(102, 202)])
      g.mag.ci <- which(mag.md %in% g.mag)
      g.mag.co <- mag.rel[, g.mag.ci]
      g.mag.yr <- yrs[g.mag.ci]
      g.mag.dr <- g.dur[g.mag.ci]
      g.mag.id <- cpt.id[g.mag.ci]
    }
  }

  if (as_list) {
    v <- list(est = y$est, cpt = cpts, slo = slo.seg, mag = mag, mag_rel = mag.rel,              # matrices
           len = len.seg, cpt_id = cpt.id, tpa = tpav, d_dur = d.dur, g_dur = g.dur,             # long vectors (lenght equal to that of the time series)
           d_max = d.mag.co, d_fst = d.fst.mg, g_max = g.mag.co, wgts = w,                       # short vectors (length equal to the number of bands)
           d_max_md = d.mag, d_fst_md = d.fst.md, g_mag_md = g.mag, d_max_yr = d.mag.yr,         # single values
           d_fst_yr = d.fst.yr, g_max_yr = g.mag.yr, d_max_dr = d.mag.dr, d_fst_dr = d.fst.dr,   # single values
           g_max_dr = g.mag.dr, d_max_id = d.mag.id, d_fst_id = d.fst.id, g_max_id = g.mag.id,   # single values
           n_gap = ngap, n_tpa = ntpa, n_iter = iter)                                            # single values
  }
  else {
    v <- c(y$est, cpts, slo.seg, mag, mag.rel,              # matrices
           len.seg, cpt.id, tpav, d.dur, g.dur,             # long vectors (lenght equal to that of the time series)
           d.mag.co, d.fst.mg, g.mag.co, w,                 # short vectors (length equal to the number of bands)
           d.mag, d.fst.md, g.mag, d.mag.yr, d.fst.yr,      # single values
           g.mag.yr, d.mag.dr, d.fst.dr, g.mag.dr,          # single values
           d.mag.id, d.fst.id, g.mag.id, ngap, ntpa, iter)  # single values
  }
}

