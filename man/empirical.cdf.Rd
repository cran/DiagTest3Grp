\name{empirical.cdf}
\alias{empirical.cdf}

\title{
 Empirical cumulative density function estimation
}
\description{
 Provide empirical cumulative density function (CDF) estimation of a sample
 vector xx at a point \eqn{c_0}.  
}
\usage{
  empirical.cdf(xx, c0)  
}

\arguments{
  \item{xx}{A numeric vector, sample data.
}
  \item{c0}{A numerical value, the value at which the CDF will be estimated.
}
}
\details{
  The CDF of xx at c0 is estimated as \eqn{\sum{I_{(xx \le c_0)}/n}} where
  n is the length of xx and \eqn{I(.)} is the indicator function giving
  1 if the criterion in the bracket is satisfied and 0 otherwise.
}
\value{
  Return a numerical value indicating the cumulative probability of xx
  at \eqn{c_0} by empirical CDF method.
}
\author{
 Jingqin Luo
}

\seealso{
  \code{\link{KernelSmoothing.cdf}}
  \code{\link{Youden3Grp}}
}

\examples{

x <- rnorm(100)

res <- empirical.cdf(xx=x,c0=0.5)
}
\keyword{univar}
\keyword{nonparametric}

