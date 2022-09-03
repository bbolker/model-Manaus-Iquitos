\name{Testpoints}
\alias{sobol.testpoints}
\alias{grid.testpoints}
\title{Sobol' low-discrepancy sequence}
\description{
  These function generate data frames containing test points.
  \code{sobol.testpoints} uses the Sobol' low-discrepancy sequence to
  disperse the grid points in an optimal manner.  \code{grid.testpoints}
  generates an ordinary regular rectangular grid.
}
\usage{
sobol.testpoints(free, fixed, n)
grid.testpoints(free, fixed, n)
}
\arguments{
  \item{n}{Number of vectors requested.}
  \item{free}{Named list of ranges of variables.}
  \item{fixed}{Named list of fixed parameters.}
}
\details{
  Specifying a variable as 'free' overrides its specification as fixed.
}
\value{
  \item{sobol.testpoints}{Returns a data frame with n 'observations' of
    the variables in \code{free} and \code{fixed}.}
  \item{grid.testpoints}{Returns a data frame with \eqn{n^d} observations,
    where d is the number of 'free' variables.}
}
\references{Press et al.}
\author{Aaron A. King (kingaa at umich dot edu)}
\examples{
free <- list(a=c(0,1),b=c(100,200));
fixed <- list(c=9,f=12);
x <- sobol.testpoints(free,fixed,100);
y <- grid.testpoints(free,fixed,10);
plot(x)
plot(y)
}
\keyword{design}
