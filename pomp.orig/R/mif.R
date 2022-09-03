mif.cooling <- function (factor, n) {
  alpha <- factor^(n-1)
  list(alpha=alpha,gamma=alpha^2)
}

powerlaw.cooling <- function (init = 1, delta = 0.1, eps = (1-delta)/2, n) {
  m <- init
  if (n <= m) {                         # linear cooling regime
    alpha <- (m-n+1)/m
    gamma <- alpha^2
  } else {                              # power-law cooling regime
    alpha <- ((n/m)^(delta+eps))/n
    gamma <- (n/m)^(delta+1)/n/n
  }
  list(alpha=alpha,gamma=gamma)
}

mif <- function (object, ... )
  stop("function 'mif' is undefined for objects of class '",class(object),"'")

mif.pomp <- function (object, Nmif = 1,
                      pars = stop("'pars' must be specified"),
                      ivps = character(0),
                      stvs = character(0),
                      start = stop("'start' must be specified"),
                      rw.sd = stop("'rw.sd' must be specified"),
                      alg.pars = stop("'alg.pars' must be specified"),
                      weighted = TRUE, tol = 1e-17, warn = TRUE, max.fail = 0, ndone = 0) {
  start.names <- names(start)
  if (is.null(start.names))
    stop("mif error: 'start' must be a named vector")
  theta <- start
  sigma <- rw.sd
  npars <- length(pars)
  if (npars == 0)
    stop("mif error: 'pars' must be a nonempty character vector")
  Nv <- length(theta)
  Nt <- ncol(object@data)
  if ((length(sigma)==1) && (sigma==0)) {
    sigma <- rep(0,Nv)
    names(sigma) <- start.names
  }
  if ((length(sigma)!=Nv) || any(!(names(sigma)%in%start.names)))
    stop("'rw.sd' must be of the same length as 'start' and their names must match")
  sigma <- sigma[start.names]
  result <- matrix(NA,nrow=Nmif+1,ncol=Nv+2)
  colnames(result) <- c('loglik','nfail',names(theta))
  result[1,] <- c(NA,NA,theta)
  if (!('cooling.factor' %in% names(alg.pars)))
    stop("you must specify the exponential cooling factor 'cooling.factor' in the 'alg.pars' argument")
  for (n in 1:Nmif) {
    cool.sched <- try(
                      mif.cooling(alg.pars$cooling.factor,ndone+n),
                      silent=T
                      )
    if (inherits(cool.sched,'try-error'))
      stop("mif error: cooling schedule error\n",cool.sched)
    sigma.n <- sigma*cool.sched$alpha
    X <- try(
             particles.pomp(object,Np=alg.pars$Np,center=theta,sd=alg.pars$CC*sigma.n),
             silent=T
             )
    if (inherits(X,'try-error'))
      stop("mif error: error in 'particles'\n",X)
    x <- try(
             pfilter.pomp.internal(
                                   object,
                                   X=X,
                                   sigma=sigma.n,
                                   pars=pars,
                                   ivps=ivps,
                                   stvs=stvs,
                                   tol=tol,
                                   warn=warn,
                                   max.fail=max.fail
                                   ),
             silent=T
             )
    if (inherits(x,'try-error'))
      stop("mif error: error in 'pfilter'\n",x)

    v <- x$pred.var[pars,,drop=F]
    
    if (weighted) {                     # MIF update rule
      v1 <- cool.sched$gamma*(1+alg.pars$CC^2)*sigma[pars]^2
      theta.hat <- cbind(theta[pars],x$filter.mean[pars,,drop=F])
      theta[pars] <- theta[pars]+apply(apply(theta.hat,1,diff)/t(v),2,sum)*v1
    } else {                            # unweighted (flat) average
      theta.hat <- x$filter.mean[pars,,drop=F]
      theta[pars] <- apply(theta.hat,1,mean)
    }
    theta[ivps] <- x$filter.mean[ivps,alg.pars$T0]
    result[n+1,-1] <- c(x$nfail,theta)
    result[n,1] <- x$loglik
  }
  new(
      "mif",
      object,
      ivps=ivps,
      pars=pars,
      stvs=stvs,
      Nmif=as.integer(Nmif),
      alg.pars=alg.pars,
      coef=theta,
      random.walk.sd=sigma,
      pred.mean=x$pred.mean,
      pred.var=x$pred.var,
      filter.mean=x$filter.mean,
      conv.rec=result,
      cond.loglik=x$cond.loglik,
      eff.sample.size=x$eff.sample.size,
      loglik=x$loglik
      )
}

mif.mif <- function (object, Nmif = object@Nmif,
                     pars = object@pars, ivps = object@ivps, stvs = object@stvs,
                     start = object@coef, rw.sd = object@random.walk.sd,
                     alg.pars = object@alg.pars,
                     weighted = TRUE, tol = 1e-17, warn = TRUE, max.fail = 0,
                     ndone = 0) {
  mif.pomp(
           object,
           Nmif=Nmif,
           pars=pars,
           ivps=ivps,
           stvs=stvs,
           start=start,
           rw.sd=rw.sd,
           alg.pars=alg.pars,
           weighted=weighted,
           tol=tol,
           warn=warn,
           max.fail=max.fail,
           ndone=ndone
           )
}

setGeneric('mif')
setMethod('mif','pomp',mif.pomp)
setMethod('mif','mif',mif.mif)

continue <- function (object, ... )
  stop("function 'continue' is undefined for objects of class '",class(object),"'")

continue.mif <- function (object, Nmif, ...) {
  ndone <- object@Nmif
  obj <- mif.mif(object,Nmif=Nmif,ndone=ndone,...)
  object@conv.rec[ndone+1,'loglik'] <- obj@conv.rec[1,'loglik']
  obj@conv.rec <- rbind(
                        object@conv.rec,
                        obj@conv.rec[-1,colnames(object@conv.rec)]
                        )
  obj@Nmif <- as.integer(ndone+Nmif)
  obj
}

setGeneric('continue')
setMethod('continue','mif',continue.mif)
