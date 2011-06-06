\name{NonParametric.VUS.var}
\alias{NonParametric.VUS.var}

\title{
  Boostrap variance of the nonparametric estimate on volumn under ROC surface (VUS)  
}
\description{
  Calculate the variance of the nonparametric estimate of volumn under
  ROC surface (VUS) through bootstrapping.
}
\usage{
  NonParametric.VUS.var(x, y, z, alpha = 0.05, NBOOT = 100)  
}

\arguments{
  \item{x}{
  A numeric vector, a diagnostic test's measurements in  the \eqn{D^-} (usually healthy
  subjects).  
}
\item{y}{
A numeric vector,  a diagnostic test's measurements in  the \eqn{D^0} (usually mildly
  diseased subjects).  
}
\item{z}{
A numeric vector,  a diagnostic test's measurements in  the \eqn{D^+} (usually severely
  diseased subjects).  
}
  \item{alpha}{
A numeric value,  confidence level for confidence interval, will result in
    (1-alpha)*100\%; Default alpha=0.05, i.e., 95\% confidence
    interval.  
}
  \item{NBOOT}{
A numeric value, total number of bootstrap samples to draw for estimate of variance of
  the VUS nonparametric estimate.
}
}
\details{
Return a numeric value, the variance estimation on the nonparametric VUS estimate obtained
 through bootstrapping.
}
\value{ Return a list of outputs, including,
  \item{variance}{A numeric value, the estimate of variance of
    nonparametric VUS through bootstrapping.} 
  \item{lower}{A numeric value, the lower limit of bootstrap basic
    quantile confidence interval for the estimate of variance of
    nonparametric VUS, i.e., (alpha/2)*100\% quantile. Under default
    alpha=0.05, lower= basic 0.25\% quantile of bootstrap variance
    estimates of the nonparametric VUS estimator.}
  \item{upper}{A numeric value, the upper limit of bootstrap basic
    quantile confidence interval for the estimate of variance of
    nonparametric VUS, i.e., (1-alpha/2)*100\% quantile. Under default
    alpha=0.05, upper= basic 0.975\% quantile of bootstrap variance
    estimates of the nonparametric VUS estimator.}
}
\references{
Bradley Efron, Robert Tibshirani, Robert J. Tibshirani (1993) An introduction
to the bootstrap \emph{CRC Press}.
}
\author{
 Jingqin Luo
}
\note{
  Bug reports, malfunctioning, or suggestions for further improvements or
 contributions can be sent to Jingqin Luo <rosy@wubios.wustl.edu>.  
}

\section{Warning}{The bootstrap variance estimate may take a minute. }
  

\seealso{
  \code{\link{VUS}}  \code{\link{Normal.VUS}}  \code{\link{Youden3Grp.Variance.Bootstrap}} 
}
\examples{

 data(AL)
 group <- AL$group
 table(group)

 ##take the negated FACTOR1 marker measurements
 factor1 <- -AL$FACTOR1

 x <- factor1[group=="D-"]
 y <- factor1[group=="D0"]
 z <- factor1[group=="D+"]

 nonpar.vus <- NonParametric.VUS.var (x=x, y=y, z=z, alpha=0.05, NBOOT=100)
}
\keyword{nonparametric}
\keyword{htest}
