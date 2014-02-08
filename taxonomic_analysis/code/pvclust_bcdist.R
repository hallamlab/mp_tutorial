pvclust <- function(data, method.hclust="average",
                    method.dist="correlation", use.cor="pairwise.complete.obs",
                    nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE)
  {
    # data: (n,p) matrix, n-samples, p-variables
    n <- nrow(data); p <- ncol(data)

    # hclust for original data
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                 "median", "centroid")
    method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
    distance <- dist.pvclust(data, method=method.dist, use.cor=use.cor)
    data.hclust <- hclust(distance, method=method.hclust)
    
    # multiscale bootstrap
    size <- floor(n*r)
    rl <- length(size)
    
    if(rl == 1) {
      if(r != 1.0)
        warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
      
      r <- list(1.0)
    }
    else
      r <- as.list(size/n)
    
    mboot <- lapply(r, boot.hclust, data=data, object.hclust=data.hclust, nboot=nboot,
                    method.dist=method.dist, use.cor=use.cor,
                    method.hclust=method.hclust, store=store, weight=weight)
    
    result <- pvclust.merge(data=data, object.hclust=data.hclust, mboot=mboot)
    
    return(result)
  }

plot.pvclust <- function(x, print.pv=TRUE, print.num=TRUE, float=0.01,
                         col.pv=c(2,3,8), cex.pv=0.8, font.pv=NULL,
                         col=NULL, cex=NULL, font=NULL, lty=NULL, lwd=NULL,
                         main=NULL, sub=NULL, xlab=NULL, ...)
{
  if(is.null(main))
    main="Cluster dendrogram with AU/BP values (%)"
  
  if(is.null(sub))
    sub=paste("Cluster method: ", x$hclust$method, sep="")
  
  if(is.null(xlab))
    xlab=paste("Distance: ", x$hclust$dist.method)
      
  plot(x$hclust, main=main, sub=sub, xlab=xlab, col=col, cex=cex,
       font=font, lty=lty, lwd=lwd, ...)

  if(print.pv)
    text(x, col=col.pv, cex=cex.pv, font=font.pv, float=float, print.num=print.num)
}

text.pvclust <- function(x, col=c(2,3,8), print.num=TRUE,  float=0.01, cex=NULL, font=NULL, ...)
{
  axes <- hc2axes(x$hclust)
  usr  <- par()$usr; wid <- usr[4] - usr[3]
  au <- as.character(round(x$edges[,"au"]*100))
  bp <- as.character(round(x$edges[,"bp"]*100))
  rn <- as.character(row.names(x$edges))
  au[length(au)] <- "au"
  bp[length(bp)] <- "bp"
  rn[length(rn)] <- "edge #"
  a <- text(x=axes[,1], y=axes[,2] + float * wid, au,
            col=col[1], pos=2, offset=.3, cex=cex, font=font)
  a <- text(x=axes[,1], y=axes[,2] + float * wid, bp,
            col=col[2], pos=4, offset=.3, cex=cex, font=font)
  if(print.num)
    a <- text(x=axes[,1], y=axes[,2], rn,
              col=col[3], pos=1, offset=.3, cex=cex, font=font)
}

print.pvclust <- function(x, which=NULL, digits=3, ...)
{
  if(is.null(which)) which <- 1:nrow(x$edges)
  cat("\n")
  cat(paste("Cluster method: ", x$hclust$method, "\n", sep=""))
  cat(paste("Distance      : ", x$hclust$dist.method, "\n\n", sep=""))
  cat("Estimates on edges:\n\n")
  print(round(x$edges[which,], digits=digits))
  cat("\n")
}

summary.pvclust <- function(object, ...){
  class(object) <- "list"
  summary(object, ...)
}

pvrect <- function(x, alpha=0.95, pv="au", type="geq", max.only=TRUE, border=2, ...)
  {
    len <- nrow(x$edges)
    member <- hc2split(x$hclust)$member
    order  <- x$hclust$order
    usr <- par("usr")
    xwd <- usr[2] - usr[1]
    ywd <- usr[4] - usr[3]
    cin <- par()$cin

    ht <- c()
    j <- 1

    if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
       stop("Invalid type argument: see help(pvrect)")
    
    for(i in (len - 1):1)
      {
        if     (pm==1) wh <- (x$edges[i,pv] >= alpha) # Greater than or EQuals
        else if(pm==2) wh <- (x$edges[i,pv] <= alpha) # Lower than or EQuals
        else if(pm==3) wh <- (x$edges[i,pv] >  alpha) # Greater Than
        else if(pm==4) wh <- (x$edges[i,pv] >  alpha) # Lower Than

        if(wh)
          {
            mi <- member[[i]]
            ma <- match(mi, order)
            
            if(max.only == FALSE || (max.only && sum(match(ma, ht, nomatch=0)) == 0))
              {
                xl <- min(ma)
                xr <- max(ma)
                yt <- x$hclust$height[i]
                yb <- usr[3]
                
                mx <- xwd / length(member) / 3
                my <- ywd / 200
                
                rect(xl - mx, yb + my, xr + mx, yt + my, border=border, shade=NULL, ...)
                
                j <- j + 1
              }
            ht <- c(ht, ma)
          }
      }
  }

msplot <- function(x, edges=NULL, ...)
  {
    if(is.null(edges)) edges <- 1:length(x$msfit)
    d   <- length(edges)

    mfrow.bak <- par()$mfrow
    on.exit(par(mfrow=mfrow.bak))

    par(mfrow=n2mfrow(d))

    for(i in edges) {
      if(i == 1 || (i %% 10 == 1 && i > 20))
        main <- paste(i, "st edge", sep="")
      else if(i == 2 || (i %% 10 == 2 && i > 20))
        main <- paste(i, "nd edge", sep="")
      else if(i == 3 || (i %% 10 == 3 && i > 20))
        main <- paste(i, "rd edge", sep="")
      else
        main <- paste(i, "th edge", sep="")

      plot(x$msfit[[i]], main=main, ...)
    }
  }

lines.pvclust <- function(x, alpha=0.95, pv="au", type="geq", col=2, lwd=2, ...)
  {
    len <- nrow(x$edges)
    member <- hc2split(x$hclust)$member
    order  <- x$hclust$order
    usr <- par("usr")
    xwd <- usr[2] - usr[1]
    ywd <- usr[4] - usr[3]
    cin <- par()$cin

    ht <- c()
    j <- 1
    
    if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
       stop("Invalid type argument: see help(lines.pvclust)")
    
    for(i in (len - 1):1)
      {
        if     (pm==1) wh <- (x$edges[i,pv] >= alpha) # Greater than or EQuals
        else if(pm==2) wh <- (x$edges[i,pv] <= alpha) # Lower than or EQuals
        else if(pm==3) wh <- (x$edges[i,pv] >  alpha) # Greater Than
        else if(pm==4) wh <- (x$edges[i,pv] >  alpha) # Lower Than

        if(wh)
          {
            mi <- member[[i]]
            ma <- match(mi, order)
            
            if(sum(match(ma, ht, nomatch=0)) == 0)
              {
                xl <- min(ma)
                xr <- max(ma)
                yt <- x$hclust$height[i]
                yb <- usr[3]
                
                mx <- xwd/length(member)/10
                
                segments(xl-mx, yb, xr+mx, yb, xpd=TRUE, col=col, lwd=lwd, ...)

                j <- j + 1
              }
            ht <- c(ht, ma)
          }
      }
  }

pvpick <- function(x, alpha=0.95, pv="au", type="geq", max.only=TRUE)
  {
    len <- nrow(x$edges)
    member <- hc2split(x$hclust)$member
    order  <- x$hclust$order
    
    ht <- c()
    a  <- list(clusters=list(), edges=c()); j <- 1

    if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
       stop("Invalid type argument: see help(pickup)")
    
    for(i in (len - 1):1)
      {
        if     (pm==1) wh <- (x$edges[i,pv] >= alpha) # Greater than or Equals
        else if(pm==2) wh <- (x$edges[i,pv] <= alpha) # Lower than or Equals
        else if(pm==3) wh <- (x$edges[i,pv] >  alpha) # Greater Than
        else if(pm==4) wh <- (x$edges[i,pv] >  alpha) # Lower Than

        if(wh)
          {
            mi <- member[[i]]
            ma <- match(mi, order)

            if(max.only == FALSE || (max.only && sum(match(ma, ht, nomatch=0)) == 0))
              {
                a$clusters[[j]] <- x$hclust$labels[mi]
                a$edges <- c(a$edges,i)
                
                j <- j + 1
              }
            ht <- c(ht, ma)
          }
      }
    
    a$edges <- a$edges[length(a$edges):1]
    a$clusters <- a$clusters[length(a$edges):1]

    return(a)
  }

parPvclust <- function(cl, data, method.hclust="average",
                       method.dist="correlation", use.cor="pairwise.complete.obs",
                       nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE,
                       weight=FALSE,
                       init.rand=TRUE, seed=NULL)
  {
    if(!(require(snow))) stop("Package snow is required for parPvclust.")
    
    if((ncl <- length(cl)) < 2 || ncl > nboot) {
      warning("Too small value for nboot: non-parallel version is executed.")
      return(pvclust(data,method.hclust,method.dist,use.cor,nboot,r,store))
    }

    if(init.rand) {
      if(is.null(seed))
        seed <- 1:length(cl)
      else if(length(seed) != length(cl))
        stop("seed and cl should have the same length.")
      
      # setting random seeds
      parLapply(cl, as.list(seed), set.seed)
    }

    # data: (n,p) matrix, n-samples, p-variables
    n <- nrow(data); p <- ncol(data)

    # hclust for original data
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
                 "median", "centroid")
    method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
    distance <- dist.pvclust(data, method=method.dist, use.cor=use.cor)
    data.hclust <- hclust(distance, method=method.hclust)
    
    # multiscale bootstrap
    size <- floor(n*r)
    rl <- length(size)
    
    if(rl == 1) {
      if(r != 1.0)
        warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
      
      r <- list(1.0)
    }
    else
      r <- as.list(size/n)

    nbl <- as.list(rep(nboot %/% ncl,times=ncl))
    
    if((rem <- nboot %% ncl) > 0)
    nbl[1:rem] <- lapply(nbl[1:rem], "+", 1)

    cat("Multiscale bootstrap... ")
    
    mlist <- parLapply(cl, nbl, pvclust.node,
                       r=r, data=data, object.hclust=data.hclust, method.dist=method.dist,
                       use.cor=use.cor, method.hclust=method.hclust,
                       store=store, weight=weight)
    cat("Done.\n")
    
    mboot <- mlist[[1]]

    for(i in 2:ncl) {
      for(j in 1:rl) {
        mboot[[j]]$edges.cnt <- mboot[[j]]$edges.cnt + mlist[[i]][[j]]$edges.cnt
        mboot[[j]]$nboot <- mboot[[j]]$nboot + mlist[[i]][[j]]$nboot
        mboot[[j]]$store <- c(mboot[[j]]$store, mlist[[i]][[j]]$store)
      }
    }

    result <- pvclust.merge( data=data, object.hclust=data.hclust, mboot=mboot)
    
    return(result)
  }

msfit <- function(bp, r, nboot) {

  if(length(bp) != length(r))
    stop("bp and r should have the same length")

  nboot <- rep(nboot, length=length(bp))

  use <- bp > 0 & bp < 1

  p <- se <- c(0,0); names(p) <- names(se) <- c("au", "bp")
  coef <- c(0,0); names(coef) <- c("v", "c")

  a <- list(p=p, se=se, coef=coef, df=0, rss=0, pchi=0); class(a) <- "msfit"

  if(sum(use) < 2) {
    # if(mean(bp) < .5) a$p[] <- c(0, 0) else a$p[] <- c(1, 1)
    if(mean(bp) < .5) a$p[] <- c(0, bp[r==1.0]) else a$p[] <- c(1, bp[r==1.0])
    return(a)
  }

  bp <- bp[use]; r <- r[use]; nboot <- nboot[use]
  zz <- -qnorm(bp)
  vv <- ((1 - bp) * bp) / (dnorm(zz)^2 * nboot)
  a$use <- use; a$r <- r; a$zz <- zz

  X   <- cbind(sqrt(r), 1/sqrt(r)); dimnames(X) <- list(NULL, c("v","c"))
  fit <- lsfit(X, zz, 1/vv, intercept=FALSE)
  a$coef <- coef <- fit$coef

  h.au <- c(1, -1); h.bp <- c(1, 1)
  
  z.au <- drop(h.au %*% coef); z.bp <- drop(h.bp %*% coef)
  a$p["au"] <- pnorm(-z.au); a$p["bp"] <- pnorm(-z.bp)
  V <- solve(crossprod(X, X/vv))
  vz.au <- drop(h.au %*% V %*% h.au); vz.bp <- drop(h.bp %*% V %*% h.bp)
  a$se["au"] <- dnorm(z.au) * sqrt(vz.au); a$se["bp"] <- dnorm(z.bp) * sqrt(vz.bp)
  a$rss <- sum(fit$residual^2/vv)
  
  if((a$df <- sum(use) - 2) > 0) {
    a$pchi <- pchisq(a$rss, lower.tail=FALSE, df=a$df)
  }
  else a$pchi <- 1.0

  return(a)
}

plot.msfit <- function(x, curve=TRUE, main=NULL, sub=NULL, xlab=NULL, ylab=NULL, ...)
{
  if(is.null(main)) main="Curve fitting for multiscale bootstrap resampling"
  if(is.null(sub))
    {
      sub  <- paste("AU = ", round(x$p["au"], digits=2),
                    ", BP = ", round(x$p["bp"], digits=2),
                    ", v = ", round(x$coef["v"], digits=2),
                    ", c = ", round(x$coef["c"], digits=2),
                    ", pchi = ", round(x$pchi, digits=2))
    }
  if(is.null(xlab)) xlab=expression(sqrt(r))
  if(is.null(ylab)) ylab=expression(z-value)
  
  a <- sqrt(x$r); b <- x$zz
  
  if(!is.null(a) && !is.null(b)) {
    plot(a, b, main=main, sub=sub, xlab=xlab, ylab=ylab, ...)
    if(curve) lines(x, ...)
  }
  else if (!is.null(a)){
    plot(0, 0, main=main, sub=sub, xlab=xlab, ylab=ylab,
         type="n", xaxt="n", yaxt="n", ...)
    a <- text(mean(a), 0, "No fitting")
  }
}

lines.msfit <- function(x, col=2, lty=1, ...) {
  v <- x$coef["v"]; c <- x$coef["c"]
  curve(v * x + c / x, add=TRUE, col=col, lty=lty)
}

summary.msfit <- function(object, digits=3, ...) {
  cat("\nResult of curve fitting for multiscale bootstrap resampling:\n\n")

  cat("Estimated p-values:\n")
  pv <- data.frame(object$p, object$se)
  names(pv) <- c("Estimate", "Std. Error"); row.names(pv) <- c("au", "bp")
  print(pv, digits=digits); cat("\n")

  cat("Estimated coefficients:\n")
  coef <- object$coef
  print(coef, digits=digits); cat("\n")

  cat(paste("Residual sum of squares: ", round(object$rss,digits=digits)),
      ",   p-value: ", round(object$pchi, digits=digits),
      " on ", object$df, " DF\n\n", sep="")
}

seplot <- function(object, type=c("au", "bp"), identify=FALSE,
                   main=NULL, xlab=NULL, ylab=NULL, ...)
  {
    if(!is.na(pm <- pmatch(type[1], c("au", "bp")))) {
      wh <- c("au", "bp")[pm]
      
      if(is.null(main))
        main <- "p-value vs standard error plot"
      if(is.null(xlab))
        xlab <- c("AU p-value", "BP value")[pm]
      if(is.null(ylab))
        ylab <- "Standard Error"
      
      plot(object$edges[,wh], object$edges[,paste("se", wh, sep=".")],
           main=main, xlab=xlab, ylab=ylab, ...)
      if(identify)
        identify(x=object$edges[,wh], y=object$edges[,paste("se", wh, sep=".")],
                 labels=row.names(object$edges))
    }
    else stop("'type' should be \"au\" or \"bp\".")
  }

# internal

hc2axes <- function(x)
{
  A <- x$merge # (n,n-1) matrix
  n <- nrow(A) + 1
  x.axis <- c()
  y.axis <- x$height
  
  x.tmp  <- rep(0,2)
  zz     <- match(1:length(x$order),x$order)

    for(i in 1:(n-1)) {
        ai <- A[i,1]

        if(ai < 0)
          x.tmp[1] <- zz[-ai]
        else
          x.tmp[1] <- x.axis[ai]
        
        ai <- A[i,2]
        
        if(ai < 0)
          x.tmp[2] <- zz[-ai]
        else
          x.tmp[2] <- x.axis[ai]

        x.axis[i] <- mean(x.tmp)
      }
  
  return(data.frame(x.axis=x.axis,y.axis=y.axis))
}

hc2split <- function(x)
  {
    A <- x$merge # (n-1,n) matrix
    n <- nrow(A) + 1
    B <- list()

    for(i in 1:(n-1)){
        ai <- A[i,1]
        
        if(ai < 0)
          B[[i]] <- -ai
        else
          B[[i]] <- B[[ai]]        
        
        ai <- A[i,2]
        
        if(ai < 0)
          B[[i]] <- sort(c(B[[i]],-ai))
        else
          B[[i]] <- sort(c(B[[i]],B[[ai]]))
      }

    CC <- matrix(rep(0,n*(n-1)),nrow=(n-1),ncol=n)
    
    for(i in 1:(n-1)){
        bi <- B[[i]]
        m <- length(bi)
        for(j in 1:m)
          CC[i,bi[j]] <- 1
      }

    split <- list(pattern=apply(CC,1,paste,collapse=""), member=B)
    
    return(split)
  }

pvclust.node <- function(x, r,...)
  {
#    require(pvclust)
    mboot.node <- lapply(r, boot.hclust, nboot=x, ...)
    return(mboot.node)
  }

boot.hclust <- function(r, data, object.hclust, method.dist, use.cor,
                        method.hclust, nboot, store, weight=F)
{ 
  n     <- nrow(data)
  size  <- round(n*r, digits=0)
  if(size == 0)
    stop("invalid scale parameter(r)")
  r <- size/n

  pattern   <- hc2split(object.hclust)$pattern
  edges.cnt <- table(factor(pattern)) - table(factor(pattern))
  st <- list()
  
  # bootstrap start
  rp <- as.character(round(r,digits=2)); if(r == 1) rp <- paste(rp,".0",sep="")
  cat(paste("Bootstrap (r = ", rp, ")... ", sep=""))
  w0 <- rep(1,n) # equal weight
  na.flag <- 0
  
  for(i in 1:nboot){
    if(weight && r>10) {  ## <- this part should be improved
      w1 <- as.vector(rmultinom(1,size,w0)) # resampled weight
      suppressWarnings(distance <- distw.pvclust(data,w1,method=method.dist,use.cor=use.cor))
    } else {
      smpl <- sample(1:n, size, replace=TRUE)
      suppressWarnings(distance  <- dist.pvclust(data[smpl,],method=method.dist,use.cor=use.cor))
    }
    if(all(is.finite(distance))) { # check if distance is valid
      x.hclust  <- hclust(distance,method=method.hclust)
      pattern.i <- hc2split(x.hclust)$pattern # split
      edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern))
    } else {
      x.hclust <- NULL
	  na.flag <- 1
    }

    if(store)
      st[[i]] <- x.hclust
  }
  cat("Done.\n")
  # bootstrap done
  
  if(na.flag == 1)
	warning(paste("inappropriate distance matrices are omitted in computation: r = ", r), call.=FALSE)

  boot <- list(edges.cnt=edges.cnt, method.dist=method.dist, use.cor=use.cor,
               method.hclust=method.hclust, nboot=nboot, size=size, r=r, store=st)
  class(boot) <- "boot.hclust"
  
  return(boot)
}

pvclust.merge <- function(data, object.hclust, mboot){
  
  pattern <- hc2split(object.hclust)$pattern
  
  r     <- unlist(lapply(mboot,"[[","r"))
  nboot <- unlist(lapply(mboot,"[[","nboot"))
  store <- lapply(mboot,"[[", "store")
  
  rl <- length(mboot)
  ne <- length(pattern)
  
  edges.bp <- edges.cnt <- data.frame(matrix(rep(0,ne*rl),nrow=ne,ncol=rl))
  row.names(edges.bp) <- pattern
  names(edges.cnt) <- paste("r", 1:rl, sep="")

  for(j in 1:rl) {
    edges.cnt[,j] <- as.vector(mboot[[j]]$edges.cnt) 
    edges.bp[,j]  <- edges.cnt[,j] / nboot[j]
  }
  
  ms.fitted <- lapply(as.list(1:ne),
                      function(x, edges.bp, r, nboot){
                        msfit(as.vector(t(edges.bp[x,])), r, nboot)},
                      edges.bp, r, nboot)
  class(ms.fitted) <- "mslist"
  
  p    <- lapply(ms.fitted,"[[","p")
  se   <- lapply(ms.fitted,"[[","se")
  coef <- lapply(ms.fitted,"[[","coef")
  
  au    <- unlist(lapply(p,"[[","au"))
  bp    <- unlist(lapply(p,"[[","bp"))
  se.au <- unlist(lapply(se,"[[","au"))
  se.bp <- unlist(lapply(se,"[[","bp"))
  v     <- unlist(lapply(coef,"[[","v"))
  cc    <- unlist(lapply(coef,"[[","c"))
  pchi  <- unlist(lapply(ms.fitted,"[[","pchi"))
  
  edges.pv <- data.frame(au=au, bp=bp, se.au=se.au, se.bp=se.bp,
                         v=v, c=cc, pchi=pchi)

  row.names(edges.pv) <- row.names(edges.cnt) <- 1:ne

  result <- list(hclust=object.hclust, edges=edges.pv, count=edges.cnt,
                 msfit=ms.fitted, nboot=nboot, r=r, store=store)

  class(result) <- "pvclust"
  return(result)
}

dist.pvclust <- function(x, method="euclidean", use.cor="pairwise.complete.obs")
{
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - cor(x, method="pearson", use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(cor(x,method="pearson",use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  else if(!is.na(pmatch(method,"uncentered"))){
    if(sum(is.na(x)) > 0){
      x <- na.omit(x)
      warning("Rows including NAs were omitted")
    }
    x  <- as.matrix(x)
    P  <- crossprod(x)
    qq <- matrix(diag(P),ncol=ncol(P))
    Q  <- sqrt(crossprod(qq))
    res <- as.dist(1 - P/Q)
    attr(res,"method") <- "uncentered"
    return(res)
  } 
  else if(!is.na(pmatch(method,"bray–curtis"))) {
	try(library("ecodist"), install.packages("ecodist"))
	bcdist(t(x)) # Niels start
  }
  else
    dist(t(x),method)
}


corw <- function(x,w,
                 use=c("all.obs","complete.obs","pairwise.complete.obs")
                 ) {
  if(is.data.frame(x)) x <- as.matrix(x)
  x <- x[w>0,,drop=F]
  w <- w[w>0]

  n <- nrow(x) # sample size
  m <- ncol(x) # number of variables
  if(missing(w)) w <- rep(1,n)
  r <- matrix(0,m,m,dimnames=list(colnames(x),colnames(x)))
  diag(r) <- 1
  use <- match.arg(use)

  pairu <- F
  if(use=="all.obs") {
    u <- rep(T,n)
  } else if(use=="complete.obs") {
    u <- apply(x,1,function(y) !any(is.na(y)))
  } else if(use=="pairwise.complete.obs") {
    pairu <- T
    ux <- is.finite(x)
  } else stop("unknown use")
  
  for(i in 1+seq(length=m-1)) {
    for(j in seq(length=i-1)) {
      if(pairu) u <- ux[,i] & ux[,j]
      wu <- w[u]; xi <- x[u,i]; xj <- x[u,j]
      ws <- sum(wu)
      if(ws > 1e-8) {
        xi <- xi - sum(wu*xi)/ws
        xj <- xj - sum(wu*xj)/ws
        vxi <- sum(wu*xi*xi)/ws
        vxj <- sum(wu*xj*xj)/ws
        if(min(vxi,vxj) > 1e-8)  {
          vxij <- sum(wu*xi*xj)/ws
          rij <- vxij/sqrt(vxi*vxj)
        } else {
          rij <- 0
        }
      } else {
        rij <- 0
      }
      r[i,j] <- r[j,i] <- rij
    }
  }
  r
}

### calculate distance by weight
distw.pvclust <- function(x,w,method="correlation", use.cor="pairwise.complete.obs")
{
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - corw(x,w, use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(corw(x,w, use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  stop("wrong method")
}
