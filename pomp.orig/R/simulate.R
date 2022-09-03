simulate.pomp.internal <- function (object, X, seed = NULL, return.as) {
  rv <- pmatch(return.as,c('pomp','obs','states','both'))[1]
  if (is.na(rv))
    stop("simulate error: unrecognized value of 'return.as': ",return.as)
  Nt <- ncol(object@data)
  Nd <- nrow(object@data)-1
  if ((!is.matrix(X)) || (is.null(rownames(X))))
    stop("simulate error: 'X' must be a matrix with rownames")
  Np <- ncol(X)
  Nv <- nrow(X)
  x <- array(NA,dim=c(Nv,Np,Nt))
  times <- c(object@t0,object@data[1,])
  if (!is.null(seed)) {
    save.seed <- .Random.seed
    set.seed(seed)
  }
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
      stop("simulate error: error in user rprocess\n",X)
    x[,,nt] <- X
  }
  dimnames(x) <- list(rownames(X),NULL,NULL)
  if (!(rv==3)) { # if we only want the states, we need not apply the measurement model
    ## simulate from the measurement model
    y <- array(NA,dim=c(Nd+1,Np,Nt))
    for (nt in 1:Nt) {
      xx <- try(
                do.call(
                        object@rmeasure,
                        c(
                          list(X=as.matrix(x[,,nt]),time=times[nt+1]),
                          object@userdata
                          )
                        ),
                silent=T
                )
      if (inherits(xx,'try-error'))
        stop("simulate error: error in user rmeasure\n",xx)
      y[1,,nt] <- times[nt+1]
      y[-1,,nt] <- xx
    }
    dimnames(y) <- list(rownames(object@data),NULL,NULL)
    y <- aperm(y,c(1,3,2))
  }
  x <- aperm(x,c(1,3,2))                # now x is Nv x Nt x Np
  if (!is.null(seed)) {
    .Random.seed <<- save.seed
  }
  switch(
         rv,
         {                              # return.as == 'pomp'
           if (Np == 1) {
             res <- as(object,'pomp')
             res@data <- y[,,1]
             res
           } else {
             lapply(
                    1:Np,
                    function (k) {
                      res <- as(object,'pomp')
                      res@data <- y[,,k] # CHECK: DO ROWNAMES FOLLOW THIS ASSIGNMENT
                      res
                    }
                    )
           }
         },
         {                              # return.as == 'obs'
           y                   # return just the array of observations
         },
         {                              # return.as == 'states'
           list(times=times[-1],states=x) # return just the array of particles and the times
         },
         {                              # return.as == 'both'
           list(times=times[-1],states=x,obs=y)  # return times, states, and observations
         }
         )
}

simulate.pomp <- function (object, nsim=1, seed = NULL, coef, return.as = c('pomp','obs','states','both')) {
  if (is.null(names(coef)) || !is.null(dim(coef)))
    stop("'coef' must be a named vector")
  X <- try(
           particles(object,Np=nsim,center=coef,sd=0),
           silent=T
           )
  if (inherits(X,'try-error'))
    stop("simulate error: error in 'particles'\n",X)
  simulate.pomp.internal(object,seed=seed,X=X,return.as)
}

simulate.mif <- function (object, nsim = 1, seed = NULL, return.as = c('pomp','obs','states','both')) {
  simulate.pomp(
                object,
                nsim=nsim,
                seed=seed,
                coef=object@coef,
                return.as=return.as
                )
}

setMethod('simulate','pomp',simulate.pomp)
setMethod('simulate','mif',simulate.mif)
