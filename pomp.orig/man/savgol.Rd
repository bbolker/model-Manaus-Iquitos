\name{savgol}
\alias{savgol}
\alias{savgol.filter}
\title{Savitsky-Golay smoothing filter}
\description{
  The Savitsky-Golay filter smooths while preserving a specified number
  of moments.  The filters provided here allow for smoothing and differentiation.
}
\usage{
savgol(x, m, left, right = left, nd = 0)
savgol.filter(m, left, right = left, nd = 0)
}
\arguments{
  \item{x}{Vector or array of data to be filtered.}
  \item{m}{Number of moments to match.}
  \item{left}{Size of window to the left.}
  \item{right}{Size of window to the right.}
  \item{nd}{Derivative to approximate.}
}
\details{
  More details here.
}
\value{
  Savitsky-Golay filter.
}
\references{Press et al.}
\author{Aaron A. King (kingaa at umich dot edu)}
\note{More notes.}
\examples{
left <- 10;
right <- 10;
x <- cumsum(matrix(data = rnorm(200), 200, 1));
y <- savgol.filter(m = 4, left = left, right = right, nd = 0);
z0 <- c(rep(NA,left),convolve(x,y,type='f'),rep(NA,right));
y <- savgol.filter(m = 4, left = left, right = right, nd = 1);
z1 <- c(rep(NA,left),convolve(x,y,type='f'),rep(NA,right));
y <- savgol.filter(m = 4, left = left, right = right, nd = 2);
z2 <- c(rep(NA,left),convolve(x,y,type='f'),rep(NA,right));
X <- ts(matrix(data = c(x,z0,z1,z2), ncol=4));
plot(X)
x <- apply(matrix(data=rnorm(200),100,2),2,cumsum)
z <- savgol(x,m=4,left,right,nd=0)
plot(c(0,100),range(x),type='n')
points(x[,1],col='red')
points(x[,2],col='blue')
lines(z[,1],col='red')
lines(z[,2],col='blue')
}
\keyword{smooth}
