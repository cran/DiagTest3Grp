\name{KernelSmoothing.cdf}
\alias{KernelSmoothing.cdf}
\title{
  Kernel smoothing cumulative density function (CDF)
}
\description{
 Calculate the kernel smoothing cumulative density function (CDF) of a
 given sample data at a user-specified value.   
}
\usage{
  KernelSmoothing.cdf(xx, c0, bw)
}

\arguments{
  \item{xx}{
  A numeric vector, sample data. 
}
  \item{c0}{
 A numeric value, the cumulative probability for which \eqn{Pr(xx \le c_0)}{P(xx<=c_0)} will be
 estimated.
}
  \item{bw}{
  A numeric value indicating the bandwidth used in the kernel smoothing density approximation.
}
}
\details{
 Kernel smoothing is a popular method to approximate a probability
 density function (PDF) or cumulative density function (CDF). Normal
 density function is conveniently used in the function as the kernel
 density and bandwidth is calculated according to the normal reference
 rule or the Sheather-Jones plug-in method in the package or can be
 specified arbitrarily by users. For a sample data xx,the cumulative CDF
 \eqn{P(xx \le c0)}{P(xx<= c0)} can be approximated by kernel smoothing method as
 \eqn{\frac{1}{n}\sum_{j=1}^{n}\Phi(\frac{xx_j-c_0}{bw})} where n is the
 sample size of xx.  
}
\value{
  Return a numeric value---the cumulative probability.
}
\references{
  
  Silverman, B.W. (1986) Density Estimation for Statistics and Data
  Analysis. Chapman & Hall.
  
  Wasserman, L. (2005) All of Statistics: A Concise Course in Statistical
  Inference. \emph{Springer}
}
\author{
 Jingqin Luo
}

\note{
 Bug reports, malfunctioning, or suggestions for further improvements or
 contributions can be sent to Jingqin Luo <rosy@wubios.wustl.edu>.
}

\seealso{
 \code{\link{BW.ref}} \code{\link{Youden3Grp}}
}
\examples{

###generate data
x <- rnorm(100,10,1.5)

##calcualte bandwidth by normal refernce rule 
bw1 <- KernelSmoothing.cdf(xx=x, c0=6, bw=0.1)
}
\keyword{smooth}
\keyword{nonparametric}

