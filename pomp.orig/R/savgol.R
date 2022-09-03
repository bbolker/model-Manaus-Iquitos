savgol.filter <- function (m, left, right = left, nd = 0) {
  if ((left < 0) || (right < 0))
    stop("bad arguments: left and right must be nonnegative")
  if (nd > m)
    stop("must have nd <= m")
  if ((left + right) < m)
    stop("must have left + right >= m")
  a <- outer(-left:right,0:m,FUN='^')
  s <- svd(a)
  s$u %*% diag(1/s$d) %*% s$v[nd+1,]
}

savgol <- function (x, m, left, right = left, nd = 0) {
  y <- savgol.filter(m,left,right,nd)
  if (is.array(x)) {
    z <- apply(x,-1,function(w)c(rep(NA,left),convolve(w,y,type='f'),rep(NA,right)))
  } else {
    z <- c(rep(NA,left),convolve(x,y,type='f'),rep(NA,right))
  }
  return(z)
}
