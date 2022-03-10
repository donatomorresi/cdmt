### misc functions
finding.cp <- function(bts.obj){

  if(length(bts.obj$sameboat)==1L){
    all.edges <- matrix(c(bts.obj$decomp.hist[1,,], bts.obj$sameboat), ncol=1) #sameboat shows pairs
  } else{
    all.edges <- rbind(bts.obj$decomp.hist[1,,], bts.obj$sameboat) #sameboat shows pairs
  }
  survived.edges <- all.edges[ , abs(bts.obj$decomp.hist[3,1,]) > sqrt(.Machine$double.eps), drop=F]

  if(length(survived.edges)==0L){
    cp <- c()
  } else if(length(survived.edges)>0L & dim(survived.edges)[2]>1L){
    i <- 1L
    cp <- c()
    while(i<dim(survived.edges)[2]){
      part.obj <- bts.obj$decomp.hist[1,c(1,2),]
      matched <- which(!is.na(match(data.frame(part.obj), data.frame(matrix(survived.edges[c(2:3),i], ncol=1)))))

      if((survived.edges[4,i]!=0L) & diff2(survived.edges[-4,i])[1]==1L & diff2(survived.edges[-4,i+1L])[1]==1L){
        cp <- c(cp, survived.edges[c(1, 3), i])
        i <- i+2L
      } else if((survived.edges[4,i]!=0L) & diff2(survived.edges[-4,i])[2]==1L & diff2(survived.edges[-4,i+1])[1]==1L){
        cp <- c(cp, survived.edges[c(1, 3), i+1L])
        i <- i+2L
      } else if((survived.edges[4,i]==0L) & diff2(survived.edges[-4,i])[1]==1L & diff2(survived.edges[-4,i])[2]!=1L){
        cp <-c(cp, survived.edges[c(1,3), i])
        i <- i+1L
      } else if((survived.edges[4,i]==0L) & diff2(survived.edges[-4,i])[1]==1L & diff2(survived.edges[-4,i])[2]==1L & length(matched)>0L){
        cp <-c(cp, survived.edges[c(1,2), i])
        i <- i+1L
      } else {
        cp <-c(cp, survived.edges[c(1,2,3), i])
        i <- i+1L
      }
    }
    cp <- unique(sort(cp))
    cp <- c(cp, bts.obj$n+1L)
  } else if(length(survived.edges)>0L & dim(survived.edges)[2]==1L & survived.edges[4,1]==0L & diff2(survived.edges[-4,1])[1]==1L & diff2(survived.edges[-4,1])[2]!=1L){
    cp <- c()
    cp <- c(cp, survived.edges[c(1,3), 1])
  } else if(length(survived.edges)>0L & dim(survived.edges)[2]==1L & survived.edges[4,1]==0L & diff2(survived.edges[-4,1])[1]==1L & diff2(survived.edges[-4,1])[2]==1L){ ### 3) (x, xx from a chunk)
    cp <- c()
    cp <- c(cp, survived.edges[c(1,2), 1])
  } else {
    cp <- c()
  }

  ### last adjustment
  if(bts.obj$n==3L & length(cp)>0L & dim(survived.edges)[2]==1L){
    cp <- cp[-1]
    cp <- c(cp, bts.obj$n)
  } else {
    cp <- cp[(cp<=bts.obj$n & cp>1L)]
  }

  ### for comparing with other methods
  if(length(cp)>0L){
    cp <- cp-1L
  }
  return(cp)
}


computeDET <- function(edges = edges, edgerow = edgerow, weights.const = weights.const, weights.lin = weights.lin, bts.coeffs = bts.coeffs){

  sub.wc <- cbind(weights.const[edges[edgerow,1]], weights.const[edges[edgerow,2]], weights.const[edges[edgerow,3]])
  sub.wl <- cbind(weights.lin[edges[edgerow,1]], weights.lin[edges[edgerow,2]], weights.lin[edges[edgerow,3]])
  sub.tc <- c(bts.coeffs[edges[edgerow,1]], bts.coeffs[edges[edgerow,2]], bts.coeffs[edges[edgerow,3]])

  sub.m <- cbind(sub.wc, sub.wl)
  detcoef <- matrix(nrow=3L, ncol=length(edgerow))

  for (q in seq_len(length(edgerow))) {
    detcoef[,q] <- filter.bts(sub.m[q,])
  }

  details <- colSums2(detcoef * matrix(sub.tc, nrow=3, byrow=T))
  dim(sub.tc) <- c(length(edgerow), 3)

  return(list(detcoef=detcoef, det=details, wc=sub.wc, wl=sub.wl, tc=sub.tc))
}


updating <- function(ee = ee, weights.const = weights.const, weights.lin = weights.lin, bts.coeffs = bts.coeffs, idx = idx){

  wc0 <- cbind(weights.const[ee[,1]], weights.const[ee[,2]], weights.const[ee[,3]])
  wl0 <- cbind(weights.lin[ee[,1]], weights.lin[ee[,2]], weights.lin[ee[,3]])
  tc0 <- cbind(bts.coeffs[ee[,1]], bts.coeffs[ee[,2]], bts.coeffs[ee[,3]])

  m <- cbind(wc0,wl0)
  h <- matrix(nrow=3L, ncol=dim(m)[1])

  for (q in seq_len(dim(h)[2])) {
    h[,q] <- filter.bts(m[q,])
  }

  M0 <- matrix(nrow=dim(h)[2], ncol=9L)

  for (k in seq_len(dim(h)[2])) {
    M0[k,] <- c(orth.matrix(h[,k]))
  }

  wc1 <- cbind(rowSums2(wc0*M0[,c(1,4,7)]), rowSums2(wc0*M0[,c(2,5,8)]), rowSums2(wc0*M0[,c(3,6,9)]))
  wl1 <- cbind(rowSums2(wl0*M0[,c(1,4,7)]), rowSums2(wl0*M0[,c(2,5,8)]), rowSums2(wl0*M0[,c(3,6,9)]))
  tc1 <- cbind(rowSums2(tc0*M0[,c(1,4,7)]), rowSums2(tc0*M0[,c(2,5,8)]), rowSums2(tc0*M0[,c(3,6,9)]))

  eating.up0 <- ee[,1]
  eating.up1 <- ee[,2]
  eaten.up <- ee[,3]
  idx <- idx[is.na(match(idx, c(eaten.up)))]

  ### 1) updating X
  bts.coeffs[eating.up0] <- c(tc1[,2])
  bts.coeffs[eating.up1] <- c(tc1[,3])
  bts.coeffs[eaten.up] <- c(tc1[,1])

  ### 2) updating weight.const
  weights.const[eaten.up] <- c(wc1[,1])
  weights.const[eating.up0] <- c(wc1[,2])
  weights.const[eating.up1] <- c(wc1[,3])

  ### 3) updating weight.lin
  weights.lin[eaten.up] <- c(wl1[,1])
  weights.lin[eating.up0] <- c(wl1[,2])
  weights.lin[eating.up1] <- c(wl1[,3])

  return(list(weights.const=weights.const, weights.lin=weights.lin, bts.coeffs=bts.coeffs, idx=idx, h=h, tc1=tc1))
}


balance.np <- function(paired=paired, ee=ee, idx, no.of.current.steps=no.of.current.steps, n=n){

  prd <- !is.na(matrix(match(ee, paired), ncol=3))
  blnc <- matrix(nrow=3, ncol=no.of.current.steps)
  firsttwo <- which(rowSums2(prd[,1:2, drop=F])==2L)
  lasttwo <- which(rowSums2(prd[,2:3, drop=F])==2L)

  if(length(c(firsttwo, lasttwo))>0L){
    nopair <- seq_len(dim(ee)[1])[-c(firsttwo, lasttwo)]
  } else{
    nopair <- seq_len(dim(ee)[1])
  }

  if(length(firsttwo)>0L){
    prtn <- ee[firsttwo,3] - ee[firsttwo,1]
    blnc[1:2, firsttwo] <- rbind(prtn/(prtn+1L), 1/(prtn+1L))
  }
  if(length(lasttwo)>0L){
    prtn <- idx[match(ee[lasttwo, 3], idx)+1L] - ee[lasttwo, 2]
    if(sum(is.na(prtn))>0L){
      prtn[is.na(prtn)] <- n - ee[is.na(prtn), 2] + 1L
    }
    blnc[1:2, lasttwo] <- rbind(1/(prtn+1L), prtn/(prtn+1L))
  }
  if(length(nopair)>0L){
    blnc[1:3, nopair] <- 1/3
  }
  return(blnc)
}


balance.p <- function(pr=pr, ee.p1=ee.p1, idx, ee.p2=ee.p2, n=n){

  blnc <- matrix(nrow=3, ncol=dim(pr)[2])

  c1 <- ee.p1[,3] - ee.p1[,1]
  c2 <- idx[match(ee.p2[, 3], idx)+1] - ee.p1[, 3]
  if(sum(is.na(c2))>0L){
    c2[is.na(c2)] <- n - ee.p1[is.na(c2), 3] + 1L
  }

  blnc[1:2,] <- matrix(c(c1/(c1+c2), c2/(c1+c2)), ncol=2, byrow=T)

  return(blnc)
}


filter.bts <- function(a) {

  a5a1_a4a2.dif <- a[5]*a[1] - a[4]*a[2]
  a2a6_a3a5.dif <- a[2]*a[6] - a[3]*a[5]
  a4a3_a1a6.dif <- a[4]*a[3] - a[1]*a[6]

  w <- -sqrt( a5a1_a4a2.dif^2 / (a2a6_a3a5.dif^2 + a4a3_a1a6.dif^2 + a5a1_a4a2.dif^2))

  u <- w * a2a6_a3a5.dif / a5a1_a4a2.dif

  v <- w * a4a3_a1a6.dif / a5a1_a4a2.dif

  df <- c(u, v, w)

  if (any(is.na(df))) {
    z <- filter.bts(a[6:1])
    df <- z[3:1]
  }
  return(df)

}


orth.matrix <- function(d) {

  M <- matrix(0L, 3, 3)

  M[1,] <- d
  M[1,] <- M[1,] / sqrt(sum(M[1,]^2))
  u <- M[1, 1]
  v <- M[1, 2]
  w <- M[1, 3]

  M[2,] <- c(1-u^2, -u*v, -u*w)
  M[3,] <- c(0, -w, v)

  M[2,] <- M[2,] / sqrt(sum(M[2,]^2))
  M[3,] <- M[3,] / sqrt(sum(M[3,]^2))

  return(M)

}


L12 <- function(l){
  x <- seq_len(l)
  n <- length(x)

  edges <- matrix(0L, n-2, 3) # to be updated for each scale j
  for(i in seq_len(n-2L)){
    edges[i,] <- c((n-i-1L):(n-i+1L))
  }

  weights.const <- rep(1L, n)
  weights.lin <- seq_len(n)
  updatedS <- diag(n)
  L1L2 <- vector('list', dim(edges)[1])

  for(st in seq_len(dim(edges)[1])) {
    ee <- matrix(edges[st,], 1, 3)

    h <- filter.bts(c(weights.const[edges[st,]], weights.lin[edges[st,]] ))
    M <- orth.matrix(h)

    tmp <- matrix(0L, nrow=3, ncol=2)
    tmp[,1] <- weights.const[ee]
    tmp[,2] <- weights.lin[ee]

    sm.det <- M %*% tmp
    updatedS[, ee[,1]] <- tcrossprod(updatedS[, edges[st,]], M[2,,drop=F])
    updatedS[, ee[,2]] <- tcrossprod(updatedS[, edges[st,]], M[3,,drop=F])
    updatedS[, ee[,3]] <- tcrossprod(updatedS[, edges[st,]], M[1,,drop=F])

    L1L2[[st]] <- updatedS[ee[,1]:n, ee[,1:2], drop=F]

    eating.up0 <- ee[,1]
    eating.up1 <- ee[,2]
    eaten.up <- ee[,3]

    weights.const[eaten.up] <- sm.det[1, 1]
    weights.const[eating.up0] <- sm.det[2, 1]
    weights.const[eating.up1] <- sm.det[3, 1]

    weights.lin[eaten.up] <- sm.det[1, 2]
    weights.lin[eating.up0] <- sm.det[2, 2]
    weights.lin[eating.up1] <- sm.det[3, 2]
  }

  return(L1L2)
}

