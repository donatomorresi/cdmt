####################################################################
### Execute CDMT at the pixel level
cdmt_int <- function(x, nb, ny, yrs, ev, dir, th_const, noise_rm) {

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

  # return results as vector
  v <- c(y$est, cpts, slo.seg, mag, mag.rel,              # matrices
         len.seg, cpt.id, tpav, d.dur, g.dur,             # long vectors (lenght equal to that of the time series)
         d.mag.co, d.fst.mg, g.mag.co, w,                 # short vectors (length equal to the number of bands)
         d.mag, d.fst.md, g.mag, d.mag.yr, d.fst.yr,      # single values
         g.mag.yr, d.mag.dr, d.fst.dr, g.mag.dr,          # single values
         d.mag.id, d.fst.id, g.mag.id, ngap, ntpa, iter)  # single values

}

