particles <- function(object, ...)
  stop("function 'particles' is undefined for objects of class '",class(object),"'")

particles.pomp <- function (object, Np, center, sd) {
  if ((length(sd)==1) && (sd == 0)) {
    sd <- rep(0,length(center))
    names(sd) <- names(center)
  }
  if (is.null(names(center)) || is.null(names(sd)))
    stop("particles error: 'center' and 'sd' must have names")
  x <- try(
           do.call(
                   object@particles,
                   c(
                     list(Np=Np,center=center,sd=sd),
                     object@userdata
                     )
                   ),
           silent=T
           )
  if (inherits(x,'try-error'))
    stop("particles error: error in user 'particles'\n",x)
  x
}

particles.mif <- function (object, Np = NULL, center = NULL, sd = NULL) {
  if (is.null(Np)) Np <- object@Np
  if (is.null(center)) center <- object@coef
  if (is.null(sd)) sd <- object@random.walk.sd
  particles.pomp(object,Np=Np,center=center,sd=sd)
}

setGeneric('particles')
setMethod('particles','pomp',particles.pomp)
setMethod('particles','mif',particles.mif)

