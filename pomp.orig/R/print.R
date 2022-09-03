print.pomp <- function (x) {
  cat(
      "\n",
      nrow(x@data)-1,
      "-variable time series of length:",
      ncol(x@data),
      "\n"
      )
  invisible(x)
}

print.mif <- function (x) {
  s <- mif.cooling(x@alg.pars$cooling.factor,n=x@Nmif)$alpha
  print.pomp(x)
  cat(
      "number of MIF iterations done:",
      x@Nmif,
      "\n",
      length(x@pars),
      " parameter(s) estimated:",
      paste(x@pars,collapse=', '),
      "\n",
      length(x@stvs),
      " state variable(s) (STV):",
      paste(x@stvs,collapse=', '),
      "\n",
      length(x@ivps),
      " initial-value parameter(s) (IVP):",
      paste(x@ivps,collapse=', '),
      "\nSigma (parameter scales):\n"
      )
  print(x@random.walk.sd,digits=3)
  cat(
      "\nnumber of particles =",
      x@alg.pars$Np,
      "\nCC =",
      x@alg.pars$CC,
      ", T0 =",
      x@alg.pars$T0,
      "\ncooling factor =",
      x@alg.pars$cooling.factor,
      "\nrandom walk intensity at last MIF iteration:\n"
      )
  print(s*x@random.walk.sd,digits=3)
  cat(
      "\nestimated parameter(s):\n"
      )
  print(coef(x)[c(x@pars,x@ivps)],digits=3)
  cat(
      "\nlog-likelihood (w/ variable parameters):",
      round(x@loglik,1),
      "\nnumber of filtering failures:",
      x@conv.rec[x@Nmif+1,'nfail'],
      "\n\n"
      )
  invisible(x)
}

setMethod('print','pomp',print.pomp)
setMethod('print','mif',print.mif)

show.pomp <- function (object) {
  print.pomp(object)
  invisible(NULL)
}

show.mif <- function (object) {
  print.mif(object)
  invisible(NULL)
}

setMethod('show','pomp',show.pomp)
setMethod('show','mif',show.mif)

