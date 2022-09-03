\name{cqlf}
\alias{cqlf}
\title{Composite quasilikelihood function.}
\description{
Compute the quasilikelihood function given a set of simulated data and an actual data set of the same length.
}
\usage{
cqlf(simdata, data, maxlag=ncol(data))
}
\arguments{
  \item{simdata}{A matrix containing the simulated data sets, each data set in its own column.}
  \item{data}{The actual data set.  This must be a vector or 1-D array.}
  \item{maxlag}{The maximum lag to be used in the computation.}
}
\details{
  Suppose the data are \eqn{\{y_t\}_{t=1}^T} and that, for given values
  of the parameters, the model predicts mean and variance \eqn{\mu_t},
  \eqn{\sigma_t^2} in \eqn{y_t}.  Let \eqn{z_t=((y_t-\mu_t)/\sigma_t)} and
  \eqn{\rho_{st}} be the expected correlation between times \eqn{s} and
  \eqn{t}.  The composite quasilikelihood \eqn{Q} is defined by \deqn{-2
    \log{Q}=\sum_{t=1}^{T}
    \sum_{s=t-m}^{t-1} \{(z_s^2+z_t^2-2 \rho_{st} z_s
    z_t)/(1-\rho_{st}^2)+\log{\sigma_s}+\log{\sigma_t}+\log{(1-\rho_{st}^2)}\}}.
}
\value{
\eqn{-2 \log{Q}}
}
\references{TBA}
\author{Aaron A. King (kingaa at umich dot edu)}
\examples{}
\keyword{models}
\keyword{ts}
