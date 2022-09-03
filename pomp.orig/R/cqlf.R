cqlf <- function(simdata, data, maxlag=ncol(data)) {
  simdata <- drop(simdata)
  if (!is.array(simdata)) stop('simdata must be an array')
  m <- dim(simdata)
  lm <- length(m)
  data <- drop(data)
  if (!is.array(data)) data <- array(data,dim=length(data))
  n <- dim(data)
  ln <- length(n)
  if (ln > 2) stop('data should be a vector or at most 2-D array')
  if (lm != ln + 1) stop('simdata should be an array with one more dimension than data')
  if (any(m[-lm] != n)) stop('the actual and simulated datasets must agree in dimension');
  if (any(is.infinite(simdata)) || any(is.infinite(data))) return(Inf)
  if (any(is.nan(simdata)) || any(is.nan(data))) return(NaN)
  if (any(is.na(simdata)) || any(is.na(data))) return(NA)
  maxlag <- pmin(maxlag,n[1])
  if (length(maxlag) < ln) maxlag <- rep(maxlag,ln)
  if (is.na(maxlag[1])) maxlag[1] <- 0;
  if (is.na(maxlag[2])) maxlag[2] <- -1;
  .C("composite_quasilikelihood",
     Q = as.double(0),
     dim = as.integer(m),
     ldim = as.integer(lm),
     x = as.double(simdata),
     y = as.double(data),
     maxlag = as.integer(maxlag),
     DUP = FALSE,
     PACKAGE = "pomp.orig"
     )$Q
}

cqlf.pomp <- function (object, Np, theta, maxlag=ncol(object@data)) {
  Nt <- ncol(object@data)
  X <- do.call(
               object@rmodel,
               c(
                 list(theta,Np),
                 object@userdata
                 )
               )
  log(cqlf(X,object@data,maxlag))
}
