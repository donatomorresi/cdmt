########################################################
### decomposition
hd.bts.dcmp <- function(x, p = .01, w) {

  d <- dim(x)[1]
  n <- dim(x)[2]
  noe <- n-2L

  weights.const <- rep(1L, n)
  weights.lin <- seq_len(n)
  idx <- seq_len(n)
  paired <- c()

  edges <- matrix(0L, noe, 3)

  edges[,1] <- seq_len(n-2L)
  edges[,2] <- 2:(n-1L)
  edges[,3] <- 3:n

  edges <- cbind(edges, 0L)

  decomp.hist <- array(0L, dim=c(4 * d, 3, n - 2L))

  bts.coeffs <- x

  steps.left <- n - 2L
  current.step <- 0L

  sameboat <- c()

  while (dim(edges)[1]) {

    max.current.steps <- ceiling(p * steps.left)
    removable.nodes <- rep(1L, max(idx))

    Dmat <- matrix(0L, nrow = d, ncol = dim(edges)[1])

    if (any(edges[, 4] > 0L)){

      pr <- matrix(which(edges[, 4] != 0L), nrow = 2)

      for (j in seq_len(d)) {
        cd1 <- computeDET(edges=edges, edgerow=pr[1,], weights.const=weights.const, weights.lin=weights.lin, bts.coeffs=bts.coeffs[j,])
        detcoef <- cd1$detcoef
        p1.details <- cd1$det

        M0 <- matrix(nrow=dim(detcoef)[2], ncol = 9L)

        for (k in seq_len(dim(detcoef)[2])) {
          M0[k,] <- c(orth.matrix(detcoef[, k]))
        }

        upd.wc <- cbind(rowSums2(cd1$wc * M0[, c(2,5,8)]), rowSums2(cd1$wc * M0[, c(3,6,9)]), weights.const[edges[pr[2,], 3]])
        upd.wl <- cbind(rowSums2(cd1$wl * M0[, c(2,5,8)]), rowSums2(cd1$wl * M0[, c(3,6,9)]), weights.lin[edges[pr[2,], 3]])
        upd.bts.coeffs <- rbind(rowSums2(cd1$tc * M0[, c(2,5,8)]), rowSums2(cd1$tc * M0[, c(3,6,9)]), bts.coeffs[j, edges[pr[2,], 3]])
        p2.m <- cbind(upd.wc, upd.wl)
        p2.d <- matrix(nrow = 3L, ncol = dim(p2.m)[1])

        for (q in seq_len(dim(p2.d)[2])) {
          p2.d[, q] <- filter.bts(p2.m[q,])
        }

        p2.details <- colSums2(p2.d * upd.bts.coeffs)

        p.detail <- rep(rowMaxs(cbind(abs(p1.details), abs(p2.details))), each = 2)

        if (length(pr) != dim(edges)[1]){

          cd3 <- computeDET(edges = edges, edgerow = seq_len(dim(edges)[1])[-c(pr)], weights.const = weights.const, weights.lin = weights.lin, bts.coeffs = bts.coeffs[j,])
          np.details <- cd3$det
          Dmat[j, ] <- c(p.detail, np.details)
        }
        else {
          Dmat[j,] <- p.detail
        }
      }

    }
    else {
      edgerow <- seq_len(dim(edges)[1])
      sub.wc <- cbind(weights.const[edges[edgerow, 1]], weights.const[edges[edgerow, 2]], weights.const[edges[edgerow, 3]])
      sub.wl <- cbind(weights.lin[edges[edgerow, 1]], weights.lin[edges[edgerow, 2]], weights.lin[edges[edgerow, 3]])
      sub.m <- cbind(sub.wc, sub.wl)
      detcoef <- matrix(nrow = 3L, ncol = dim(sub.m)[1])

      for (q in seq_len(dim(detcoef)[2])) {
        detcoef[, q] <- filter.bts(sub.m[q,])
      }

      bts.coeffs.nrow <- dim(bts.coeffs)[1]
      detcoef.ncol <- dim(detcoef)[2]

      Dmat <- bts.coeffs[, edges[, 1], drop = FALSE] * matrix(detcoef[1,], bts.coeffs.nrow, detcoef.ncol, byrow=TRUE) +
              bts.coeffs[, edges[, 2], drop = FALSE] * matrix(detcoef[2,], bts.coeffs.nrow, detcoef.ncol, byrow=TRUE) +
              bts.coeffs[, edges[, 3], drop = FALSE] * matrix(detcoef[3,], bts.coeffs.nrow, detcoef.ncol, byrow=TRUE)

    }

    if (is.null(w)) {
      Dmat2 <- abs(Dmat)
    }
    else {
      Dmat2 <- abs(Dmat) * w
    }

    colmaxD <- colMaxs(Dmat2) + colMeans2(Dmat2)

    ord.det <- order(colmaxD)

    cand <- rbind(ord.det, edges[ord.det,4])

    eitr <- 1L
    tei <- 1L
    if(cand[2,1] > 0L){
      removable.nodes[edges[ord.det[1:2],1]] <- removable.nodes[edges[ord.det[1:2],2]] <- removable.nodes[edges[ord.det[1:2],3]] <- 0L
      tei <- tei + 1L
      eitr <- c(eitr, tei)
    } else{
      removable.nodes[edges[ord.det[1],1]] <- removable.nodes[edges[ord.det[1],2]] <- removable.nodes[edges[ord.det[1],3]] <- 0L
    }

    while ((length(eitr) < max.current.steps) & (tei < noe)) {

      tei <- tei + 1L

      if (cand[2, tei] > 0L) {
        if (sum(removable.nodes[edges[ord.det[tei:(tei+1)],1]])==2L & sum(removable.nodes[edges[ord.det[tei:(tei+1)],2]])==2L & sum(removable.nodes[edges[ord.det[tei:(tei+1)],3]])==2L){
          removable.nodes[edges[ord.det[tei:(tei+1L)],1]] <- removable.nodes[edges[ord.det[tei:(tei+1L)],2]] <- removable.nodes[edges[ord.det[tei:(tei+1L)],3]] <- 0L
          eitr <- c(eitr, tei, tei+1L)
          tei <- tei + 1L

        }
      }
      else {
        if(removable.nodes[edges[ord.det[tei],1]] & removable.nodes[edges[ord.det[tei],2]] & removable.nodes[edges[ord.det[tei],3]]){
          eitr <- c(eitr, tei)
          removable.nodes[edges[ord.det[tei],1]] <- removable.nodes[edges[ord.det[tei],2]] <- removable.nodes[edges[ord.det[tei],3]] <- 0L
        }
      }
    }

    details.min.ind <- ord.det[eitr]

    no.of.current.steps <- length(eitr)

    ee <- matrix(edges[details.min.ind,], no.of.current.steps, 4)
    sameboat <- c(sameboat, c(ee[,4]))
    idx0 <- idx

    if (sum(ee[,4]>0L)==0L){
      ee <- ee[, -4, drop=F]

      for (j in seq_len(d)) {
        udt <- updating(ee=ee, weights.const=weights.const, weights.lin=weights.lin, bts.coeffs=bts.coeffs[j,], idx=idx0)
        bts.coeffs[j,] <- udt$bts.coeffs

        decomp.hist[(j-1)*4+1,,(current.step+1):(current.step+no.of.current.steps)] <- t(ee)
        decomp.hist[(j-1)*4+2,,(current.step+1):(current.step+no.of.current.steps)] <- udt$h
        decomp.hist[(j-1)*4+3,,(current.step+1):(current.step+no.of.current.steps)] <- t(udt$tc1)
        decomp.hist[(j)*4,,(current.step+1):(current.step+no.of.current.steps)] <- balance.np(paired=paired, ee=ee, idx=idx0, no.of.current.steps=no.of.current.steps, n=n)
      }
      weights.const <- udt$weights.const
      weights.lin <- udt$weights.lin
      idx <- udt$idx

    }
    else {

      pr <- matrix(which(ee[,4]!=0L), nrow=2)
      ee[pr[2,], 1:2] <- ee[pr[1,], 1:2]
      ee <- ee[, -4, drop=F]
      ee.p1 <- ee[pr[1,],,drop=F]
      ee.p2 <- ee[pr[2,],,drop=F]
      ee.np <- ee[-c(pr),,drop=F]

      for (j in seq_len(d)) {

        udt <- updating(ee=ee.p1, weights.const=weights.const, weights.lin=weights.lin, bts.coeffs=bts.coeffs[j,], idx=idx0)
        bts.coeffs[j,] <- udt$bts.coeffs
        wg.c <- udt$weights.const
        wg.l <- udt$weights.lin
        idx <- udt$idx

        decomp.hist[(j-1)*4+1,,(current.step+1):(current.step+dim(ee.p1)[1])] <- t(ee.p1)
        decomp.hist[(j-1)*4+2,,(current.step+1):(current.step+dim(ee.p1)[1])] <- udt$h
        decomp.hist[(j-1)*4+3,,(current.step+1):(current.step+dim(ee.p1)[1])] <- t(udt$tc1)
        decomp.hist[(j-1)*4+4,,(current.step+1):(current.step+dim(ee.p1)[1])] <- balance.p(pr=pr, ee.p1=ee.p1, idx=idx0, ee.p2=ee.p2, n=n)

        udt <- updating(ee=ee.p2, weights.const=wg.c, weights.lin=wg.l, bts.coeffs=bts.coeffs[j,], idx=idx)
        bts.coeffs[j,] <- udt$bts.coeffs
        wg.c <- udt$weights.const
        wg.l <- udt$weights.lin
        idx <- udt$idx

        decomp.hist[(j-1)*4+1,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- t(ee.p2)
        decomp.hist[(j-1)*4+2,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- udt$h
        decomp.hist[(j-1)*4+3,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- t(udt$tc1)
        decomp.hist[(j-1)*4+4,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- balance.p(pr=pr, ee.p1=ee.p1, idx=idx, ee.p2=ee.p2, n=n)

        if (length(pr)!=dim(ee)[1]){

          udt <- updating(ee=ee.np, weights.const=wg.c, weights.lin=wg.l, bts.coeffs=bts.coeffs[j,], idx=idx)
          bts.coeffs[j,] <- udt$bts.coeffs

          decomp.hist[(j-1)*4+1,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- t(ee.np)
          decomp.hist[(j-1)*4+2,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- udt$h
          decomp.hist[(j-1)*4+3,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- t(udt$tc1)
          decomp.hist[(j-1)*4+4,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- balance.np(paired=paired, ee=ee.np, idx=idx,
          no.of.current.steps=no.of.current.steps-length(pr), n=n)
        }
      }

      weights.const <- udt$weights.const
      weights.lin <- udt$weights.lin
      idx <- udt$idx

    }

    paired <- sort(unique(c(paired, c(ee[,1:2]))))

    if(length(which(!is.na(match(paired, c(ee[,3])))))>0L){
      paired <- sort(paired[is.na(match(paired, c(ee[,3])))])
    }

    edges <- matrix(0L, length(idx)-2L, 3)
    edges[,1] <- idx[seq_len(length(idx)-2L)]
    edges[,2] <- idx[2:(length(idx)-1L)]
    edges[,3] <- idx[3:(length(idx))]

    matchpair <- matrix(match(edges, paired), ncol=3)
    matchpair[!is.na(matchpair)] <- 1L

    rs <- which(rowSums2(matchpair, na.rm=T)==3L)

    if(length(rs)>0L){
      edges <- as.matrix(rbind(edges[rs,], edges[-rs,]))
      matchpair <- rbind(matchpair[rs,], matchpair[-rs,])
      edges <- as.matrix(cbind(edges, c(rep(seq_len(length(rs)/2), each=2), rep(0L, dim(edges)[1]-length(rs)))))
    } else {
      edges <- as.matrix(cbind(edges, rep(0L, dim(edges)[1])))
    }

    removed <- c(which(rowSums2(matchpair, na.rm=T)==1L), which(matchpair[,1]==1L & is.na(matchpair[,2]) & matchpair[,3]==1L))
    if(length(removed)>0L){
      edges <- edges[-removed,,drop=F]
    }

    noe <- dim(edges)[1]
    steps.left <- steps.left - no.of.current.steps
    current.step <- current.step + no.of.current.steps

    if(noe==1L & dim(edges)[1]==1L){
      edges[,4] <- 0L
    }

  }

  return(list(n = n, sameboat=sameboat, decomp.hist=decomp.hist, bts.coeffs=bts.coeffs))

}

#########################################################
### Thresholding
hd.bts.dns <- function(bts.obj, lambda, bal = 0, w, agg) {

  if(dim(bts.obj$decomp.hist)[1]>4){
    d <- dim(bts.obj$decomp.hist)[1]/4    # n time series
    detail.all <- bts.obj$decomp.hist[c(4*seq_len(d)-1L),1,]

    if (agg == 'max') {
      details <- colMaxs(abs(detail.all))
    }
    else {
      details <- colWeightedMeans(abs(detail.all), w)
    }

  } else {
    d <- 1L
    detail.all <- bts.obj$decomp.hist[3, 1, , drop=FALSE]
    dim(detail.all) <- c(d, prod(dim(detail.all)[2:3]))
    details <- abs(bts.obj$decomp.hist[3, 1,])
  }

  sameboat <- bts.obj$sameboat
  n <- bts.obj$n

  protected <- rep(0L, n)

  for (i in seq_len(n-2L)) {

    if (!protected[bts.obj$decomp.hist[1,1,i]] & !protected[bts.obj$decomp.hist[1,2,i]] &
        !protected[bts.obj$decomp.hist[1,3,i]])

      bts.obj$decomp.hist[c(4*seq_len(d)-1L),1,i] <- bts.obj$decomp.hist[c(4*seq_len(d)-1L),1,i] * (
        (details[i] > lambda) &
          (bts.obj$decomp.hist[4,1,i] > bal) &
          (bts.obj$decomp.hist[4,2,i] > bal)
      )

    if (abs(bts.obj$decomp.hist[3,1,i]) > 0L) protected[bts.obj$decomp.hist[1,1,i]] <- protected[bts.obj$decomp.hist[1,2,i]] <- 1L

  }

  paired <- matrix(which(bts.obj$sameboat!=0L), nrow=2)

  if(length(paired)>0L){
    for(i in seq_len(dim(paired)[2])){

      overzero <- is.element(paired[,i], which(abs(bts.obj$decomp.hist[3,1,])>0L))
      zero <- is.element(paired[,i], which(abs(bts.obj$decomp.hist[3,1,])==0L))
      if(sum(overzero)==1L & sum(zero)==1L){
        bts.obj$decomp.hist[c(4*seq_len(d)-1L),1,paired[overzero==F,i]] <- detail.all[,paired[overzero==F,i]]
      }
    }
  }

  return(bts.obj)

}

#########################################################
### inverse transformation
hd.bts.inv <- function(bts.obj) {

  n <- bts.obj$n
  d <- dim(bts.obj$decomp.hist)[1]/4

  for (i in (n-2L):1L) {

    M.inv <- t(orth.matrix(bts.obj$decomp.hist[2,,i]))

    ind <- bts.obj$decomp.hist[1,,i]

    bts.obj$decomp.hist[c(4*seq_len(d)-1L),2,i] <- bts.obj$bts.coeffs[,ind[1]]
    bts.obj$decomp.hist[c(4*seq_len(d)-1L),3,i] <- bts.obj$bts.coeffs[,ind[2]]

    tmp <- bts.obj$decomp.hist[c(4*seq_len(d)-1L),,i]

    if(d==1L){
      rcstr.tmp <- M.inv %*% matrix(tmp, ncol=1)
    } else{
      #rcstr.tmp <- M.inv %*% t(tmp)
      rcstr.tmp <- tcrossprod(M.inv, tmp)
    }

    bts.obj$bts.coeffs[,ind[1]] <- rcstr.tmp[1,]
    bts.obj$bts.coeffs[,ind[2]] <- rcstr.tmp[2,]
    bts.obj$bts.coeffs[,ind[3]] <- rcstr.tmp[3,]

  }

  return(bts.obj)

}

#########################################################
### post processing - stage 1
hd.bts.pp1 <- function(bts.obj, lambda, w) {

  wc <- rep(1L, bts.obj$n)
  wl <- seq_len(bts.obj$n)
  pp1fit <- bts.obj$bts.coeffs
  inicp <- finding.cp(bts.obj)

  if(length(inicp) > 0L){
    chp <- c(1L, inicp+1L, bts.obj$n+1L)
    pqr <- cbind(chp[seq_len(length(chp)-2)], chp[2:(length(chp)-1L)]-1L, chp[3:length(chp)]-1L)
    d.pqr <- rowDiffs(pqr)

    details <- detail.1 <- detail.2 <- matrix(nrow=dim(pp1fit)[1], ncol=dim(d.pqr)[1])

    while(length(chp)>2L){
      ######################################
      ### 1) (x, x)
      c1 <- which(d.pqr[,1]==0L & d.pqr[,2]==1L)
      if(length(c1)>0L){
        for(i in seq_along(c1)){
          p <- pqr[c1[i],1]
          q <- pqr[c1[i],2]
          r <- pqr[c1[i],3]
          details[,c1[i]] <- abs((pp1fit[,p:q, drop=F]-pp1fit[,(q+1L):r, drop=F])/sqrt(2))
        }
      }
      ######################################
      ### 23) (x, x, x)
      c23 <- which((d.pqr[,1]==0L & d.pqr[,2]==2L) | (d.pqr[,1]==1L & d.pqr[,2]==1L))
      if(length(c23)>0L){
        for(i in seq_along(c23)){
          p <- pqr[c23[i],1]
          r <- pqr[c23[i],3]
          h <- filter.bts(c(wc[p:r], wl[p:r]))
          details[,c23[i]] <- abs(pp1fit[,p:r, drop=F]%*%matrix(h, ncol=1))
        }
      }
      ######################################
      ### 4) (xx, xx)
      c4 <- which(d.pqr[,1]==1L & d.pqr[,2]==2L)
      if(length(c4)>0L){
        for(i in seq_along(c4)){
          p1 <- pqr[c4[i],1]
          q1 <- pqr[c4[i],2]
          r1 <- q1+1L
          h1 <- filter.bts(c(wc[p1:r1], wl[p1:r1]))
          detail.1[,c4[i]] <- abs(pp1fit[,p1:r1, drop=F]%*%matrix(h1, ncol=1))

          ### updating
          M <- orth.matrix(h1)

          r2 <- pqr[c4[i],3]
          h2 <- filter.bts(c((M%*%matrix(wc[p1:r1], ncol=1))[2:3,1], wc[r2], (M%*%matrix(wl[p1:r1], ncol=1))[2:3,1], wl[r2]))
          detail.2[,c4[i]] <- abs(cbind(t(tcrossprod(M, pp1fit[,p1:r1, drop=F]))[,2:3,drop=F], pp1fit[, r2, drop=F])%*%matrix(h2, ncol=1))

          details[,c4[i]] <- rowMaxs(cbind(abs(detail.1[,c4[i]]), abs(detail.2[,c4[i]])))
        }
      }
      ######################################
      ### 56) (x, xx from chunk) or (xx from chunk, x)
      c5 <- which(d.pqr[,1]==0L & d.pqr[,2]>2L)
      if(length(c5)>0L){
        for(i in seq_along(c5)){
          p <- pqr[c5[i],1]
          q <- pqr[c5[i],2]
          r <- pqr[c5[i],3]

          l12 <- L12(r-q)[[r-q-2L]]
          wcL <- matrix(wc[(q+1L):r], nrow=1)%*%l12
          wlL <- matrix(wl[(q+1L):r], nrow=1)%*%l12
          xL <- pp1fit[,(q+1L):r, drop=F]%*%l12

          h <- filter.bts(c(wc[p:q], wcL , wl[p:q], wlL))
          details[,c5[i]] <- abs(cbind(pp1fit[,p:q, drop=F], xL)%*%matrix(h, ncol=1))
        }
      }

      c6 <- which(d.pqr[,1]>1L & d.pqr[,2]==1L)
      if(length(c6)>0L){
        for(i in seq_along(c6)){
          p <- pqr[c6[i],1]
          q <- pqr[c6[i],2]
          r <- pqr[c6[i],3]

          l12 <- L12(q-p+1L)[[q-p+1L-2L]]
          wcL <- matrix(wc[p:q], nrow=1)%*%l12
          wlL <- matrix(wl[p:q], nrow=1)%*%l12
          xL <- pp1fit[,p:q, drop=F]%*%l12

          h <- filter.bts(c(wcL, wc[(q+1):r], wlL, wl[(q+1L):r]))
          details[,c6[i]] <- abs(cbind(xL, pp1fit[,(q+1L):r, drop=F])%*%matrix(h, ncol=1))
        }
      }
      ######################################
      ### 78) (xx, xx from chunk) or (xx from chunk, xx)
      c7 <- which(d.pqr[,1]==1L & d.pqr[,2]>2L)
      if(length(c7)>0L){
        for(i in seq_along(c7)){
          p <- pqr[c7[i],1]
          q <- pqr[c7[i],2]
          r <- pqr[c7[i],3]

          l12 <- L12(r-q)[[r-q-2L]]
          wcL <- matrix(wc[(q+1L):r], nrow=1)%*%l12
          wlL <- matrix(wl[(q+1L):r], nrow=1)%*%l12
          xL <- pp1fit[,(q+1L):r, drop=F]%*%l12

          new.wc <- c(wc[p:q], wcL)
          new.wl <- c(wl[p:q], wlL)
          new.x <- cbind(pp1fit[,p:q, drop=F], xL)

          ### first edge
          h1 <- filter.bts(c(new.wc[1:3], new.wl[1:3]))
          detail.1[,c7[i]] <- abs(new.x[,1:3,drop=F]%*%matrix(h1, ncol=1))

          ### updating
          M <- orth.matrix(h1)

          ### second edge
          h2 <- filter.bts(c( (M%*%matrix(new.wc[1:3], ncol=1))[2:3,1], new.wc[4], (M%*%matrix(new.wl[1:3], ncol=1))[2:3,1], new.wl[4]))
          detail.2[,c7[i]] <- abs(cbind(t(tcrossprod(M, new.x[,1:3,drop=F]))[,2:3,drop=F], new.x[,4])%*%matrix(h2, ncol=1))

          details[,c7[i]] <- rowMaxs(cbind(abs(detail.1[,c7[i]]), abs(detail.2[,c7[i]])))
        }
      }

      c8 <- which(d.pqr[,1]>1L & d.pqr[,2]==2L)
      if(length(c8)>0L){
        for(i in seq_along(c8)){
          p <- pqr[c8[i],1]
          q <- pqr[c8[i],2]
          r <- pqr[c8[i],3]

          l12 <- L12(q-p+1L)[[q-p+1L-2L]]
          wcL <- matrix(wc[p:q], nrow=1)%*%l12
          wlL <- matrix(wl[p:q], nrow=1)%*%l12
          xL <- pp1fit[,p:q, drop=F]%*%l12

          new.wc <- c(wcL, wc[(q+1L):r])
          new.wl <- c(wlL, wl[(q+1L):r])
          new.x <- cbind(xL, pp1fit[,(q+1L):r, drop=F])

          ### first edge
          h1 <- filter.bts(c(new.wc[1:3], new.wl[1:3]))
          detail.1[,c8[i]] <- abs(new.x[,1:3,drop=F]%*%matrix(h1, ncol=1))

          ### updating
          M <- orth.matrix(h1)

          ### second edge
          h2 <- filter.bts(c( (M%*%matrix(new.wc[1:3], ncol=1))[2:3,1], new.wc[4], (M%*%matrix(new.wl[1:3], ncol=1))[2:3,1], new.wl[4]))
          detail.2[,c8[i]] <- abs(cbind(t(tcrossprod(M, new.x[,1:3,drop=F]))[,2:3,drop=F], new.x[,4])%*%matrix(h2, ncol=1))

          details[,c8[i]] <- rowMaxs(cbind(abs(detail.1[,c8[i]]), abs(detail.2[,c8[i]])))
        }
      }
      ######################################
      ### 9) the others
      all.idx <- seq_len(dim(d.pqr)[1])
      c9 <- all.idx[!all.idx %in% unique(c(c1, c23, c4, c5, c6, c7, c8))]

      if(length(c9) > 0L){
        for(i in seq_along(c9)){
          p <- pqr[c9[i],1]
          q <- pqr[c9[i],2]
          r <- pqr[c9[i],3]

          l12 <- L12(q-p+1L)[[q-p+1L-2L]]
          wcL1 <- matrix(wc[p:q], nrow=1)%*%l12
          wlL1 <- matrix(wl[p:q], nrow=1)%*%l12
          xL1 <- pp1fit[,p:q, drop=F]%*%l12

          l12 <- L12(r-q)[[r-q-2L]]
          wcL2 <- matrix(wc[(q+1L):r], nrow=1)%*%l12
          wlL2 <- matrix(wl[(q+1L):r], nrow=1)%*%l12
          xL2 <- pp1fit[,(q+1L):r, drop=F]%*%l12

          new.wc <- c(wcL1, wcL2)
          new.wl <- c(wlL1, wlL2)
          new.x <- cbind(xL1, xL2)

          ### first edge
          h1 <- filter.bts(c(new.wc[1:3], new.wl[1:3]))
          detail.1[,c9[i]] <- abs(new.x[,1:3,drop=F]%*%matrix(h1, ncol=1))

          ### updating
          M <- orth.matrix(h1)

          ### second edge
          h2 <- filter.bts(c( (M%*%matrix(new.wc[1:3], ncol=1))[2:3,1], new.wc[4], (M%*%matrix(new.wl[1:3], ncol=1))[2:3,1], new.wl[4]))
          detail.2[,c9[i]] <- abs(cbind(t(tcrossprod(M, new.x[,1:3,drop=F]))[,2:3,drop=F], new.x[,4])%*%matrix(h2, ncol=1))

          details[,c9[i]] <- rowMaxs(cbind(abs(detail.1[,c9[i]]), abs(detail.2[,c9[i]])))

        }
      }

      ##### Remove the change point having a smallest |detail|
      maxdet <- colMaxs(details)

      maxdetmin.i <- which.min(maxdet)

      if(maxdet[maxdetmin.i] < lambda) {
        ### update
        chp <- chp[-c(maxdetmin.i+1L)]

        for(k in seq_len(length(chp)-1L)){
          domain <- c(chp[k]:(chp[k+1L]-1L))
          X <- matrix(c(rep(1L, length(domain)), domain), ncol=2)

          for (q in seq_len(dim(pp1fit)[1])) {
            pp1fit.v <- pp1fit[q, domain]
            pp1fit[q, domain] <- pp1fit.v - .lm.fit(X, pp1fit.v)$residuals
          }
        }
      } else {
        break
      }

      ### update
      pqr <- cbind(chp[seq_len(length(chp)-2L)], chp[2:(length(chp)-1L)]-1L, chp[3:length(chp)]-1L)
      d.pqr <- rowDiffs(pqr)
      details <- detail.1 <- detail.2 <- matrix(nrow=dim(pp1fit)[1], ncol=dim(d.pqr)[1])
    }

    ### final change points
    bts.obj$chp <- chp[-c(1L, length(chp))]
    bts.obj$details <- details
    bts.obj$bts.coeffs <- pp1fit

  } else {

    details <- matrix(0L, nrow=dim(pp1fit)[1], ncol=2)

    bts.obj$chp <- inicp
    bts.obj$details <- details
    bts.obj$bts.coeffs <- pp1fit

  }

  return(bts.obj)

}

#########################################################
### all four steps (from post processing)
hd.bts.cpt <- function(x, sd, th.const, weight, p = .01, bal = 0) {

  m <- dim(x)[1]
  n <- dim(x)[2]

  x <- x/sd

  dcmp <- hd.bts.dcmp(x, p=p, w=weight)

  lambda <- th.const * sqrt(2 * log(m*n))

  agg <- c('max', 'mean')

  dns.l <- inv.l <- pp1.l <- list(length(agg))

  for (i in seq_along(agg)) {

    dns <- hd.bts.dns(dcmp, lambda=lambda, bal=bal, w=weight, agg=agg[i])
    inv <- hd.bts.inv(dns)
    pp1 <- hd.bts.pp1(bts.obj=inv, lambda=lambda, w=weight)

    dns.l[[i]] <- dns
    inv.l[[i]] <- inv
    pp1.l[[i]] <- pp1
  }

  # Compute BIC
  bic.m <- matrix(NA, m, 2)

  for (i in seq_along(agg)) {

    bts.coeffs <- pp1.l[[i]]$bts.coeffs
    resid <- x - bts.coeffs

    # log likelihood
    ll <- numeric(m)
    w0 <- rep(1,n)
    for (k in seq_len(m)) {
      ll[k] <- 0.5 * (sum(log(w0)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w0 * resid[k,]^2))))
    }

    # degrees of freedom using the penalty term proposed by Hall et al. (2013)
    ncpt <- length(pp1.l[[i]]$chp)
    df <- 2 * (ncpt + 1L) + 3*ncpt
    bic.m[,i] <- -2 * ll + log(n) * df
  }

  agg.i <- which.min(colWeightedMeans(bic.m, weight))
  pp1 <- pp1.l[[agg.i]]

  cpt <- pp1$chp

  ##cptind <- ifelse(pp1$details > lambda, 1, 0)
  cptind <- matrix(0L, dim(pp1$details)[1], dim(pp1$details)[2])
  cptind[pp1$details > lambda] <- 1L

  no.of.cpt <- length(cpt)

  est <- pp1$bts.coeffs*sd

  list(est=est, no.of.cpt=no.of.cpt, cpt=cpt, cptind=cptind)

}
