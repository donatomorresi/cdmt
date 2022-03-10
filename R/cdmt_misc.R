##################################################################
### Band-to-band correlation-based weights
wr2 <- function(x) {

  n <- dim(x)[1]
  t <- dim(x)[2]

  if (n == 1) {
    w <- 1L
    return(w)
  }
  else {

    # compute the coefficient of determination between each time series
    r2m <- matrix(nrow=n, ncol=n)

    for (i in seq_len(n)) {

      # extract the i-th time series
      y <- x[i,]

      # total sum of squares
      tss <- sum((y - mean(y))^2)

      for (j in seq_len(n)) {

        X <- matrix(c(rep(1L, t), x[j,]), ncol=2)
        res <- .lm.fit(X, y)$residuals

        # residual sum of squares
        rss <- sum(res^2)

        # coefficient of determination
        r2 <- 1 - (rss/tss)

        r2m[i,j] <- r2
      }
    }

    diag(r2m) <- 0

    # Compute the mean R2
    w <- rowSums2(r2m)/(n-1)
  }
}

##################################################################
### Estimate standard deviation using MAD
sd_est <- function(x) {

  n <- dim(x)[2] - 2L

  mad.bias <- 1 + n.times.eBias.of.mad[n] / n

  sd <- (rowMads(rowDiffs(x, differences=2)) / mad.bias) / sqrt(6)
}

##################################################################
### Detect candidate changepoints that may contain impulsive noise
cpt_cnd <- function(y, x, sd) {

  pa0 <- pa1 <- cpt.pa <- integer(0)

  if (y$no.of.cpt > 0L) {

    # exclude CPTs at extremes of the TS (not PAs)
    cpt <- y$cpt[y$cpt != dim(y$est)[2]]

    # when sequential CPTs occur, the second CPT is a PA
    seq.cpt <- diff2(cpt)

    if (any(seq.cpt == 1L)) {

      pa0 <- cpt[which(seq.cpt == 1L)]
      cpt0 <- cpt[!(cpt %in% pa0)]
      cpt.pa <- cpt0[(cpt0 - 1L) %in% pa0]
      cpt <- cpt[!(cpt %in% c(pa0, cpt.pa))]
    }

    # check for the presence of changepoints that are potential point anomalies
    if (length(cpt) > 0L) {

      # compute time offsets
      cpt.off <- unique(c(cpt, cpt - 1L))

      # compute sequential magnitudes and check for anomalies
      cnd <- logical(length(cpt.off))

      for (i in seq_along(cpt.off)) {
        cpti <- cpt.off[i]
        cnd1 <- abs(x[,cpti-1L] - x[,cpti]) > 4*sd
        cnd2 <- abs(x[,cpti-1L] - x[,cpti+1L]) < sd
        cnd[i] <- any(cnd1 & cnd2)
      }

      if (any(cnd)) {
        pa1 <- cpt.off[cnd]
      }
    }
  }

  pa <- c(pa0, pa1)

  return(pa)
}

##################################################################
### Process candidate changepoints iteratively to find noisy observations
proc_cpt <- function(x, sd, y, cpt, th_const, w, dir) {

  tpa <- integer(length(cpt))

  # discriminate between point anomalies caused by disturbances and those caused by noise
  if (length(cpt) > 0L) {

    for (i in seq_along(cpt)) {

      cpti <- cpt[i]

      # check the direction of change
      if (all(sign(x[, cpti - 1] - x[, cpti]) == -dir)) {

        tpa[i] <- cpti

      }
      else {

        # perform hits without candidate time point
        xp <- x[, -cpti, drop=F]

        # perform hits to find change points
        y <- hd.bts.cpt(xp, sd, th_const, w)

        if (any(y$cpt %in% c(cpti, cpti - 1L, cpti + 1L))) {

          if (sum(sign(y$est[, cpti-1] - y$est[, cpti]) == dir) > dim(x)[1]/2) {
            next
          }
          else {
            tpa[i] <- cpti
          }
        }
        else {
          tpa[i] <- cpti
        }
      }
    }

    if (any(tpa != 0L)) {
      tpa <- tpa[tpa != 0L]
    }
  }
  return(tpa)
}
