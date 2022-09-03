bspline <- function (x, i, degree, knots)
  .C(
     name="bspline",
     y=as.double(rep(0,length(x))),
     x=as.double(x),
     nx=as.integer(length(x)),
     i=as.integer(i),
     p=as.integer(degree),
     knots=as.double(knots),
     nknots=as.integer(length(knots)),
     DUP=FALSE,
     PACKAGE="pomp.orig"
     )$y

bspline.basis <- function (x, degree = 3, knots) {
  nbasis <- length(knots)-degree-1
  y <- matrix(0,length(x),nbasis)
  for (i in 1:nbasis) 
    y[,i] <- bspline(x,i-1,degree,knots)
  y
}

periodic.bspline.basis <- function (x, nbasis, degree = 3, period = 1) {
  if (any(x < 0) || any(x > period))
    stop("cannot evaluate the basis outside the fundamental domain")
  if (nbasis < degree)
    stop("must have nbasis < degree")
  dx <- period/nbasis
  knots <- seq(-degree*dx,period+degree*dx,by=dx)
  y <- bspline.basis(x,degree,knots)
  if (degree > 0)
    y[,1:degree] <- y[,1:degree]+y[,-(1:nbasis)]
  shift <- floor((degree-1)/2)
  y[,c((shift+1):nbasis,1:shift)]
}

