\name{pomp}
\alias{pomp}
\title{State-space models.}
\description{
  Create a new "pomp" object.
}
\usage{
  pomp(data, t0, particles, rprocess, dmeasure, rmeasure, ...)
}
\arguments{
  \item{data}{Matrix of time-series data.  The first row must correspond
    to the times at which the measurements were made.  The second and
    subsequent rows are the observations.}
  \item{t0}{The initial time for the state variables.}
  \item{particles}{Function of prototype \code{particles(Np,center,sd,...)}
    which sets up the initial particle matrix by drawing a sample of
    size \code{Np} from the initial particle distribution centered at
    \code{center} and of width \code{sd}.}
  \item{rprocess}{Function of prototype \code{rprocess(X,t1,t2,...)}
    which simulates from the unobserved process \eqn{X(t2)|X(t1)=x}.}
  \item{dmeasure}{Function of prototype \code{dmeasure(X,Y,...)} which
    gives the likelihood \eqn{f(Y=y|X=x)}.}
  \item{rmeasure}{Function of prototype \code{rmeasure(X,time,...)}
    which simulates from the observation process \eqn{Y|X=x}.}
  \item{\dots}{Any additional arguments are stored in a slot 'userdata'
    and are passed as arguments to each of the above functions.  Using
    this mechanism, the user can store additional information necessary
    for the definition of the model.}
}
\details{
  \strong{It is the user's responsibility to ensure that the
    \code{rprocess}, \code{dmeasure}, \code{rmeasure}, and
    \code{particles} elements satisfy the following conditions:}
  
  \code{rprocess} has at least the following arguments: \code{X},
  \code{t1}, \code{t2}, and \code{\dots}.  Additional arguments may be
  passed.  These will be filled with the corresponding elements added
  by the user to the 'pomp' object (see \code{pomp}).

  \code{rprocess} returns a rank-2 array (matrix) with rownames:
  each column represents a distinct particle, each row a distinct parameter
  or state variable.  The rownames are used by the algorithms (see
  \code{mif}, \code{pfilter}) to sort state variables from ordinary
  parameters from initial-value parameters.

  \code{dmeasure} has at least the arguments \code{X}, \code{Y}, and
  \code{\dots}. \code{X} will be the matrix of particles (see above).
  \code{Y} will contain a single column of the data matrix (see
  above). In particular, the first element of \code{Y} will be the time
  of the observation.  \code{dmeasure} should return a numeric vector
  containing the likelihoods of the data given each particle.
    
  \code{rmeasure} has at least the arguments \code{X}, \code{time}, and
  \code{\dots}. Additional arguments may be passed.  \code{X} must be a
  matrix with dimension \code{(nv,np)}.  Here, \code{nv} is the
  number of variables (state-variables, parameters, and IVPs) defining
  a particle and \code{np} is the number of particles.  \code{time}
  will be the time corresponding to the measurement.  \code{rmeasure}
  must return a matrix of dimension \code{(nd-1,np)}, where \code{nd} is
  the number of rows in the data matrix.  That is, if \code{y} is the
  array returned by \code{rmeasure}, then for each \code{k},
  \code{y[,,k]} is comparable to the data matrix without its first
  (time) row.

  \code{particles} has at least the following arguments: \code{Np},
  \code{center}, \code{sd}, and \code{\dots}.  Additional arguments may be
  passed.  These will be filled with the corresponding elements added
  by the user to the "pomp" object.

  \code{particles} returns a rank-2 array (matrix) with rownames:
  each column represents a distinct particle, each row a distinct
  parameter, IVP, or state variable.  The rownames are used by the
  algorithms (see \code{mif}, \code{pfilter}) to sort state variables
  from ordinary parameters from initial-value parameters.  (The user may
  choose to use the rownames in the \code{rprocess}, \code{rmeasure},
  and \code{dmeasure} functions as well.)

  The center of the particle distribution returned by \code{particles}
  should be \code{center}.  The width of the particle distribution
  should vary monotonically with \code{sd}.  In particular, when
  \code{sd=0}, \code{particles} should return a matrix consisting of copies
  of \code{center}.
}
\value{An object of class "pomp".}
\references{}
\author{Aaron A. King (kingaa at umich dot edu)}
\section{Warning}{If the user does not ensure that the model satisfies
  the above restrictions, the algorithms will probably produce
  nonsense!}
\seealso{\code{\link{pomp-class}}.  For examples, see the vignette.}
\keyword{models}
\keyword{ts}
