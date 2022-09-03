systematic.resample <- function (w) {        # systematic sampling
  Np <- length(w)
  perm <- rep(0,Np)
  .C("systematic_resampling",
     n=as.integer(Np),
     w=as.double(w),
     perm=as.integer(perm),
     DUP=TRUE,
     PACKAGE="pomp.orig"
     )$perm
}

pfilter <- function(object, ...)
  stop("function 'pfilter' is undefined for objects of class '",class(object),"'")

pfilter.pomp.internal <- function (object, X, sigma, pars, ivps, stvs, tol, warn, max.fail) {

  Nt <- ncol(object@data)
  if (is.null(dim(X))) stop("'X' must be a matrix")
  Np <- ncol(X)
  Nv <- nrow(X)

  if (is.null(rownames(X)))
    stop("pfilter error: 'X' must have rownames")
  if (any(!(c(stvs,pars,ivps)%in%rownames(X))))
    stop("pfilter error: the rownames of 'X' must include all of the parameters, ivps, and states of the pomp object 'object'")
  if (any(!(pars%in%names(sigma))))
    stop("pfilter error: the names of 'sigma' must include all of the parameters of the pomp object 'object'")

  npars <- length(pars)
  nivps <- length(ivps)
  nstvs <- length(stvs)
  req.means.vars <- (npars+nivps+nstvs > 0)
  
  loglik <- rep(NA,Nt)
  eff.sample.size <- rep(NA,Nt)

  if (req.means.vars) {
    pred.m <- matrix(
                     data=0,
                     nrow=nstvs,
                     ncol=Nt,
                     dimnames=list(stvs,NULL)
                     )
    pred.v <- matrix(
                     data=0,
                     nrow=nstvs+npars+nivps,
                     ncol=Nt,
                     dimnames=list(c(stvs,pars,ivps),NULL)
                     )
    filt.m <- matrix(
                     data=0,
                     nrow=nstvs+npars+nivps,
                     ncol=Nt,
                     dimnames=list(c(stvs,pars,ivps),NULL)
                     )
  }

  nfail <- 0

  times <- c(object@t0,object@data[1,])
  
  for (nt in 1:Nt) {

    ## advance the state variables according to the process model
    X <- try(
             do.call(
                     object@rprocess,
                     c(
                       list(X=X,t1=times[nt],t2=times[nt+1]),
                       object@userdata
                       )
                     ),
             silent=T
             )
    if (inherits(X,'try-error'))
      stop("pfilter error: error in user 'rprocess'\n",X)
        
    ## prediction means and variances
    if (req.means.vars) {
      xx <- try(
                apply(X[stvs,,drop=F],1,mean),
                silent=T
                )
      if (inherits(xx,'try-error')) {
        stop("pfilter error: error in prediction mean computation\n",xx)
      } else {
        pred.m[,nt] <- xx
      }
      xx <- try(
                apply(X[c(stvs,ivps),,drop=F],1,var),
                silent=T
                )
      if (inherits(xx,'try-error')) {
        stop("pfilter error: error in STV-IVP prediction variance computation\n",xx)
      } else {
        pred.v[c(stvs,ivps),nt] <- xx
      }
      xx <- try(
                apply(X[pars,,drop=F],1,var)+sigma[pars]^2,
                silent=T
                )
      if (inherits(xx,'try-error')) {
        stop("pfilter error: error in PAR prediction variance computation\n",xx)
      } else {
        pred.v[pars,nt] <- xx
      }
    }

    ## determine the weights
    weights <- try(
                   do.call(
                           object@dmeasure,
                           c(
                             list(X=X,Y=object@data[,nt]),
                             object@userdata
                             )
                           ),
                   silent=T
                   )
    if (inherits(weights,'try-error'))
      stop("pfilter error: error in user dmeasure\n",weights)

    ## test for failure to filter
    failures <- weights < tol
    all.fail <- all(failures)
    if (all.fail) {                # all particles are lost
      if (warn)
        message("filtering failure at time t = ",times[nt+1])
      nfail <- nfail+1
      if (nfail > max.fail)
        stop('pfilter error: too many filtering failures')
      loglik[nt] <- log(tol)          # worst log-likelihood
      weights <- rep(1/Np,Np)
      eff.sample.size[nt] <- 0
    } else {                            # not all particles are lost
      loglik[nt] <- log(mean(weights))  # compute log-likelihood
      weights[failures] <- 0
      weights <- weights/sum(weights)
      eff.sample.size[nt] <- 1/(weights%*%weights) # effective sample-size
    }

    ## compute filtering means
    if (req.means.vars) {
      filt.m[,nt] <- X[c(stvs,pars,ivps),] %*% weights
    }

    ## Matrix with samples (columns) from filtering distribution theta.t | Y.t
    if (!all.fail)
      X <- X[,systematic.resample(weights),drop=F]

    ## random walk for parameters
    if (npars > 0) {
      X[pars,] <- matrix(
                         rnorm(
                               n=Np*npars,
                               mean=X[pars,],
                               sd=sigma[pars]
                               ),
                         npars,
                         Np
                         )
    }
  }

  if (req.means.vars) {
    retval <- list(
                   loglik=sum(loglik),
                   nfail=nfail,
                   eff.sample.size=eff.sample.size,
                   cond.loglik=loglik,
                   pred.mean=pred.m,
                   pred.var=pred.v,
                   filter.mean=filt.m
                   )
  } else {
    retval <- list(
                   loglik=sum(loglik),
                   nfail=nfail,
                   eff.sample.size=eff.sample.size,
                   cond.loglik=loglik
                   )
  }

  retval
}

pfilter.pomp <- function (object,
                          Np = stop("'Np' must be specified"),
                          coef = stop("'coef' must be specified"),
                          tol = 1e-17,
                          warn = TRUE,
                          max.fail = 0
                          ) {
  sigma <- rep(0,length(coef))         # no random walk!
  names(sigma) <- names(coef)
  X <- try(
           particles.pomp(object,Np=Np,center=coef,sd=sigma),
           silent=T
           )
  if (inherits(X,'try-error'))
    stop("pfilter error: error in 'particles'\n",X)
  pfilter.pomp.internal(
                        object,
                        X=X,
                        sigma=sigma,
                        pars=character(0),
                        ivps=character(0),
                        stvs=character(0),
                        tol=tol,
                        warn=warn,
                        max.fail=max.fail
                        )
}

pfilter.mif <- function (object,
                         Np = object@alg.pars$Np,
                         coef = object@coef,
                         ...
                         ) {
  pfilter.pomp(object,Np=Np,coef=coef,...)
}

setGeneric('pfilter')
setMethod('pfilter','pomp',pfilter.pomp)
setMethod('pfilter','mif',pfilter.mif)
