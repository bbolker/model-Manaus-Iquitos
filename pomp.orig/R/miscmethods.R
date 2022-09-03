setMethod('logLik','mif',function(object)object@loglik)

setMethod('coef','mif',function(object)object@coef)

setMethod(
          'coef<-',
          'mif',
          function (object, names = NULL, ..., value) {
            if (is.null(names)) names <- names(value)
            if (is.null(names)) {
              if (length(value)!=length(object@coef))
                stop(
                     "in 'coef<-': ",
                     "number of items to replace is not equal to replacement length",
                     call.=F
                     )
              names <- names(value) <- names(object@coef)
            }
            if (is.null(names(value))) {
              if (length(value)!=length(names))
                stop(
                     "in 'coef<-': ",
                     "number of items to replace is not equal to number of names given",
                     call.=F
                     )
              names(value) <- names
            }
            incl <- names%in%names(value)
            if (!all(incl))
              stop(
                   "in 'coef<-': ",
                   "names '",paste(names[!incl],collapse=','),"' are not present in 'value'",
                   call.=F
                   )
            incl <- names%in%names(object@coef)
            if (!all(incl))
              warning(
                      "in 'coef<-': ",
                      "values '",
                      paste(names[!incl],collapse=','),
                      "' correspond to no element in 'coef(object)' and hence are not replaced",
                      call.=F
                      )
            object@coef[names[incl]] <- value[names[incl]]
            object
          }
          )

setMethod('time','pomp',function (x, ...) x@data[1,])

pred.var <- function (object, ...)
  stop("function 'pred.var' is undefined for objects of class '",class(object),"'")

pred.var.mif <- function (object, pars = NULL) {
  if (is.null(pars)) pars <- rownames(object@pred.var)
  object@pred.var[pars,]
}

setGeneric('pred.var')  
setMethod('pred.var','mif',pred.var.mif)

pred.mean <- function (object, ...)
  stop("function 'pred.mean' is undefined for objects of class '",class(object),"'")

pred.mean.mif <- function (object, pars = NULL) {
  if (is.null(pars)) pars <- rownames(object@pred.mean)
  object@pred.mean[pars,]
}

setGeneric('pred.mean')  
setMethod('pred.mean','mif',pred.mean.mif)

filter.mean <- function (object, ...)
  stop("function 'filter.mean' is undefined for objects of class '",class(object),"'")

filter.mean.mif <- function (object, pars = NULL, ...) {
  if (is.null(pars)) pars <- rownames(object@filter.mean)
  object@filter.mean[pars,]
}

setGeneric('filter.mean')  
setMethod('filter.mean','mif',filter.mean.mif)

conv.rec <- function (object, ...)
  stop("function 'conv.rec' is undefined for objects of class '",class(object),"'")

conv.rec.mif <- function (object, pars = NULL) {
  if (is.null(pars)) pars <- colnames(object@conv.rec)
  object@conv.rec[,pars]
}

setGeneric('conv.rec')  
setMethod('conv.rec','mif',conv.rec.mif)



