predvarplot.mif <- function (object, pars = object@pars, type = 'l', mean = FALSE, ...) {
  npv <- object@pred.var[pars,]/(object@random.walk.sd[pars]^2)
  nm <- rownames(object@data)
  if (is.null(nm))
    nm <- c('t',paste('observable',1:(nrow(object@data)-1)))
  if (!is.null(dim(npv))) npv <- t(npv)
  if (mean && !is.null(dim(npv)))
    npv <- apply(npv,1,mean)
  if (!is.null(dim(npv))) {
    matplot(object@data[1,],npv,type=type,ylab='prediction variance',xlab=nm[1],...)
    legend(x='topright',legend=pars,col=1:length(pars),lty=1:length(pars),bty='n')
  } else {
    plot(object@data[1,],npv,type=type,ylab='prediction variance',xlab=nm[1],...)
  }
}

plot.pomp <- function (x, y, type = 'l', ...) {
  if (dev.interactive() && (nrow(x@data) > 2)) {
    oldpar <- par(ask=T)
    on.exit(par(oldpar))
  }
  nm <- rownames(x@data)
  if (is.null(nm))
    nm <- c('t',paste('observable',1:(nrow(x@data)-1)))
  for (p in 1:(nrow(x@data)-1))
    plot(x@data[1,],x@data[p+1,],xlab=nm[1],ylab=nm[p+1],type=type,...)
}

plot.mif <- function (x, y, type = 'l', ...) {
  if (dev.interactive()) {
    oldpar <- par(ask=T)
    on.exit(par(oldpar))
  }
  plot(as(x,'pomp'),y,type=type,main='data',...)
  for (p in c('loglik',x@pars,x@ivps))
    plot(x@conv.rec[,p],main='convergence plot',ylab=p,xlab='iteration',type=type,...)
  nm <- rownames(x@data)
  if (is.null(nm))
    nm <- c('t',paste('observable',1:(nrow(x@data)-1)))
  plot(x@data[1,],x@eff.sample.size,main='effective sample size',ylab='',xlab=nm[1],type=type,...)
  plot(x@data[1,],x@cond.loglik,main='conditional log likelihood',ylab=expression(log(L)),xlab=nm[1],type=type,...)
  for (p in c(x@pars,x@stvs))
    plot(x@data[1,],x@filter.mean[p,],main='filter mean',ylab=p,xlab=nm[1],type=type,...)
}

setMethod('plot','pomp',plot.pomp)
setMethod('plot','mif',plot.mif)

predvarplot <- function (object, ...)
  stop("function 'predvarplot' is undefined for objects of class '",class(object),"'")

setGeneric('predvarplot')  
setMethod('predvarplot','mif',predvarplot.mif)

