\name{pomp-class}
\docType{class}
\alias{pomp-class}
\alias{plot,pomp-method}
\alias{plot-pomp}
\alias{print,pomp-method}
\alias{print-pomp}
\alias{simulate,pomp-method}
\alias{simulate-pomp}
\alias{mif,pomp-method}
\alias{mif-pomp}
\alias{pfilter,pomp-method}
\alias{pfilter-pomp}
\alias{show,pomp-method}
\alias{show-pomp}
\alias{time,pomp-method}
\alias{time-pomp}
\alias{particles,pomp-method}
\alias{particles-pomp}

\title{State-space model with data}
\description{The class \code{pomp} encodes a state-space model with data.}
\section{Objects from the Class}{
  Objects should be created by calls of the function
  \code{pomp}.  See the documentation for \code{\link{pomp}} for
  usage instructions and important warnings.
}
\section{Slots}{
  \describe{
    \item{\code{data}}{Matrix of time-series data.  The first row must
      correspond to the times at which the measurements were made.  The
      second and subsequent rows are the observations.}
    \item{\code{t0}:}{The initial time for the state variables.}
    \item{\code{particles}:}{Function of prototype
      \code{particles(Np,center,sd,...)} which returns a sample of size
      \code{Np} from the initial particle distribution centered at
      \code{center} and of width \code{sd}.}
    \item{\code{rprocess}:}{Function of prototype
      \code{rprocess(X,t1,t2,...)} which simulates from the unobserved
      process \eqn{X(t2)|X(t1)=x}.}
    \item{\code{dmeasure}:}{Function of prototype
      \code{dmeasure(X,Y,...)} which gives the likelihood
      \eqn{f(Y=y|X=x)}.}
    \item{\code{rmeasure}:}{Function of prototype
      \code{rmeasure(X,time,...)} which simulates from the observation
      process \eqn{Y(time)|X(time)=x}.}
    \item{\code{userdata}:}{A "list" object, containing additional data
      that will be passed to the functions above.}
  }
}
\section{Methods}{
  \describe{
    \item{mif}{Runs the MIF algorithm on the "pomp" object.  The
      call has many arguments that are described below under
      \cite{Running the MIF algorithm}.}
    \item{pfilter}{Runs a particle filter.  See the documentation below
      under \cite{Running the Particle Filter}.}
    \item{particles}{Sets up the particle matrix by making one or more
      draws from the initial particle distribution. See the
      documentation below under \cite{Initial Distribution of the
	Particles}.}
    \item{simulate}{Simulates a data series.  See the documentation
      below under \cite{Running Simulations}.}
    \item{plot}{\code{plot(object,...)} plots the data.  Additional
      arguments are passed to the low-level plotting routine.}
    \item{print}{\code{print(object,...)} prints the "pomp" object in a
      nice way.}
    \item{show}{\code{show(object)} displays the "pomp" object.}
    \item{time}{\code{time(object)} returns the vector of observation
      times.}
  }    
}
\section{Running the MIF algorithm}{
  The call sequence for the \code{mif} method is

  \code{mif(object, Nmif, pars, ivps = character(0),
    stvs = character(0), start, rw.sd, alg.pars, weighted = TRUE,
    tol = 1e-17, warn = TRUE, max.fail = 0)}

  Description of arguments:
  \describe{
    \item{object}{An object of class \code{pomp}}
    \item{Nmif}{The number of MIF iterations}
    \item{pars}{Character vector of names of ordinary parameters.}
    \item{ivps}{Character vector of names of initial-value parameters.}
    \item{stvs}{Character vector of names of state variables.}
    \item{start}{The initial guess of the parameters}
    \item{rw.sd}{The intensity of the random walk}
    \item{alg.pars}{A named list of algorithm parameters.  At the
      moment, this consists of
      \item{Np}{the number of particles to use in filtering}
      \item{CC}{the scaling coefficient relating the width of the
	initial particle distribution to \code{rw.sd}}
      \item{T0}{the timepoint for fixed-lag smoothing of initial-value
	parameters (IVPs)}
      \item{cooling.factor}{the exponential cooling factor, \eqn{a},
	where \eqn{0<a<1}}
    }
    \item{weighted}{Should a weighted average be used?}
    \item{tol}{Particles with log likelihood below \code{tol} are
      considered to be "lost".  A filtering failure occurs when, at some
      time point, all particles are lost.}
    \item{warn}{Should a warning be generated when a filtering failure
      occurs?}
    \item{max.fail}{Maximum number of filtering failures permitted.  If
      the number of failures exceeds this number, execution will
      terminate with an error.}
    }
}
\section{Running the Particle Filter}{
  One can estimate the likelihood of a constant-parameter state-space
  model using a particle filter (AKA sequential Monte Carlo) method.
  The call sequence is

  \code{pfilter(object, Np, coef, tol = 1e-17, warn = TRUE, max.fail = 0)}

  The arguments are as follows
  \describe{
    \item{object}{the "pomp" object}
    \item{Np}{the number of particles to use}
    \item{coef}{the parameters at which the log-likelihood is desired}
    \item{tol}{Particles with log likelihood below \code{tol} are
      considered to be "lost".  A filtering failure occurs when, at some
      time point, all particles are lost.}
    \item{warn}{Should the user be warned about filtering failures?}
    \item{max.fail}{The maximum number of filtering failures allowed.
      If the number of filtering failures exceeds this number,
      execution will terminate with an error.}
  }

  \code{pfilter} returns a list with the following elements
  \describe{
    \item{loglik}{The estimated log-likelihood}
    \item{nfail}{The number of filtering failures encountered}
    \item{eff.sample.size}{A vector containing the effective number of
      particles at each time point.}
    \item{cond.loglik}{A vector containing the conditional
      log-likelihoods at each time point.}
  }
}
\section{Running Simulations}{
  The call sequence for the \code{simulate} method is

  \code{simulate(object,nsim=1,seed=NULL,coef,return.as=c('pomp','obs','states','both'))}
    
  The arguments are as follows
  \describe{
    \item{object}{an object of class "pomp"}
    \item{nsim}{the number of simulations to perform}
    \item{seed}{the random seed to use}
    \item{coef}{the parameters at which to simulate the model}
    \item{return.as}{the type of object to return.  If
      \code{return.as=='pomp'}, an "pomp" object will be returned.  It
      will be identical to \code{object}, but its \code{data} slot will
      have  been replaced by the simulated data.  If \code{nsim>1}, then
      a list of "pomp" objects will be returned.  If
      \code{return.as=='obs'}, a rank-3 array with third dimension
      \code{nsim} is returned.  If \code{y} is the array returned, then
      for each \code{k}, \code{y[,,k]} is a simulated data set.  If
      \code{return.as=='states'}, a rank-3 array with the simulated
      trajectories of the \code{nsim} particles is returned.  The
      dimensions of the array is (\code{nv x nsim x nt} where \code{nv} is
      the number of variables defining a particle and \code{nt} is the
      number of times.  If \code{return.as=='both'}, a named list is
      returned.  The "obs" element, as before, contains a matrix of
      simulated observations.
    }
  }
}
\section{Initial Distribution of the Particles}{
  The \code{particles} method sets up the particle matrix by drawing
  from the initial particle distribution.  The call sequence is

  \code{particles(object, Np, center, sd)}

  The arguments are as follows
  \describe{
    \item{object}{The "pomp" object}
    \item{Np}{The number of particles, i.e., number of draws}
    \item{center}{The central value of the initial particle
      distribution}
    \item{sd}{The width of the distribution}
  }

  \code{particles} returns the particle matrix of \code{Np} columns.
}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\seealso{\code{\link{pomp}}, \code{\link{mif-class}}}
\keyword{models}
\keyword{ts}
