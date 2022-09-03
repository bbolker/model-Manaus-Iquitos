sobol <- function (dim, n) {
  x <- .C("dsobol",
          data = as.double(matrix(0,dim,n)),
          dim = as.integer(dim),
          n = as.integer(n),
          DUP = FALSE,
          PACKAGE = "pomp.orig"
          )$data;
  dim(x) <- c(dim,n);
  return(t(x));
}

sobol.testpoints <- function (free, fixed = NULL, n) {
  if (!is.list(free) || is.null(names(free))) 
    stop("'free' must be a named list");
  if (!is.null(fixed) && (!is.list(fixed) || is.null(names(fixed))))
    stop("'fixed' must be a named list");
  nm <- names(free);
  varnames <- sort(unique(c(names(fixed),nm)));
  y <- as.list(rep(NA,length(varnames)));
  names(y) <- varnames;
  fixed[nm] <- NULL;           # free overrides fixed
  y[names(fixed)] <- fixed;
  oo <- match(nm, names(y));
  if (any(is.na(oo))) 
    stop("some named arguments in 'free' are not OK");
  bds <- sapply(free, eval.parent)[,order(oo)];
  nm <- nm[order(oo)];
  d <- length(free);
  x <- sobol(d,n);
  if (d > 1) {
    for (k in 1:d) {
      y[[nm[k]]] <- bds[1,k] + (bds[2,k]-bds[1,k])*x[,k];
    }
  } else {
    y[[nm]] <- bds[1] + (bds[2]-bds[1])*x;
  }
  return(as.data.frame(y));
}

grid.testpoints <- function (free, fixed = NULL, n) {
  if (!is.list(free) || is.null(names(free))) 
    stop("'free' must be a named list");
  if (!is.null(fixed) && (!is.list(fixed) || is.null(names(fixed))))
    stop("'fixed' must be a named list");
  nm <- names(free);
  varnames <- sort(unique(c(names(fixed),nm)));
  y <- as.list(rep(NA,length(varnames)));
  names(y) <- varnames;
  fixed[nm] <- NULL;           # free overrides fixed
  y[names(fixed)] <- fixed;
  oo <- match(nm, names(y));
  if (any(is.na(oo))) 
    stop("some named arguments in 'free' are not OK");
  bds <- sapply(free, eval.parent)[,order(oo)];
  nm <- nm[order(oo)];
  d <- length(free);
  cyc <- 1:d;
  if (d > 1) {
    for (k in 1:d) {
      x <- seq(bds[1,k],bds[2,k],length=n);
      y[[nm[k]]] <- as.numeric(aperm(array(x,dim=rep(n,d)),cyc));
      cyc <- c(cyc[d], cyc[1:(d-1)]);
    }
  } else {
    y[[nm]] <- seq(bds[1],bds[2],length=n);
  }
  return(as.data.frame(y));
}
