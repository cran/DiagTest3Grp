\name{BW.ref}
\alias{BW.ref}

\title{
  Calculate the bandwidth (normal reference rule and Sheather-Jones
  plug-in rule) in kernel smoothing for density estimation
}
\description{
  
 Calculate the bandwidth in kernel smoothing for density estimation
 by use of normal reference rule and by  Sheather-Jones plug-in method.  
}
\usage{
BW.ref(x, method = "KS-SJ")
}

\arguments{
  \item{x}{
    a numeric vector, e.g., a diagnostic test's measurments.    }
  \item{method}{
    specify which method to be used for bandwidth calculation for kernel
    smoothing density estimation in the extended Youden index
    estimation. The method choices include: \dQuote{KS}--uses the simple
    normal reference rule; \dQuote{KS-SJ}--uses the Sheather-Jones plug-in
    method which used the 2-stage solve-the-equation rules calling the R
    package \emph{KernSmooth} (see references).
  }
}
\details{
 Kernel smoothing estimation of density function requires specification
 of kernel function and bandwidth. Normal density is conveniently used
 as the kernel. Bandwidth selection is important in density estimation. 

 The normal reference rule calculates the bandwidth as
 \eqn{1.06min(sd,\frac{iqr}{1.34})n^{-0.2}}{1.06*min(sd,iqr/1.34)n^{-0.2}} with
 sd and iqr indicating the standard deviation and inter-quartile range.
 The Sheather and Jones direct plug-in bandwidth estimator (JRSSB, 1991)
 is obtained by using the R package \emph{KernSmooth}.
 
}
\value{
  Return a numeric value, the calculated bandwidth for kernel smoothing
  estimate on sample data x.
}
\references{

  Sheather, S.J. and Jones, M.C. (1991) A reliable data-based bandwidth
  selection method for kernel density estimation. \emph{Journal of Royal
  Statistical Society (Serie B)} \bold{53} 683--690
  
  Sheather, S.J. (1992) The performance of six popular bandwidth selection
  methods on some real datasets. \emph{Computational Statistics}
  \bold{7} 225--250.
  
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
 \code{\link{Youden3Grp}} \code{\link{KernelSmoothing.cdf}}
}
\examples{

###generate data
x <- rnorm(100, 10, 1.5)

##calcualte bandwidth by normal refernce rule 
bw1 <- BW.ref(x=x, method="KS")

##calcualte bandwidth by Sheather-Jones plug-in method 
bw2 <- BW.ref(x=x, method="KS-SJ")
}
\keyword{univar}
\keyword{htest}
