\name{CI.normal}
\alias{CI.normal}
\title{
  Provide the lower and upper limit of a normal confidence interval
  given an estimate and its associated variance
}
\description{
 Provide the lower and upper limit of a 100*(1-alpha)\% normal
 confidence interval (CI) given an estimate and its associated variance
}
\usage{
  CI.normal(est, var0, alpha = 0.05)  
}

\arguments{
  \item{est}{
  A numeric value, a estimate on which the CI to be calculated.
}
  \item{var0}{
 A numeric value, the variance of the estimate.
}
\item{alpha}{
A numeric value, for a (1-alpha)*100\% confidence interval. Default alpha=0.05, i.e., 95\% confidence interval.
}
 }
\details{
The lower and upper limit of the confidence interval for an estimate \eqn{\hat{\theta}}{\hat{theta}}is calculated as
\eqn{\hat{\theta} -/+ Z_{\frac{\alpha}{2}}*\sqrt(var0)}{\hat{\theta} +/-
  Z_{\frac{\alpha}{2}}*\sqrt(var0)}. \eqn{Z_{\frac{\alpha}{2}}} is the normal \eqn{\alpha/2} quantile.  
}

\value{ Returns a list including,
  \item{est}{The estimate provided by user.  }
  \item{var0}{The variance of the estimate provided by user.  }
  \item{lower}{the lower limit of a 100*(1-alpha)\% normal confidence
    interval.   } 
  \item{upper}{the upper limit of a 100*(1-alpha)\% normal confidence
    interval.  } 
}
\author{
  Jingqin Luo
}

\seealso{
  \code{\link{Youden3Grp}}
  \code{\link{Youden3Grp.Variance.Bootstrap}}
}

\examples{

est <- 5

var0 <- 1.5

res <- CI.normal(est,var0,alpha=0.05)
}

\keyword{htest }

