\name{FisherZ}
\alias{FisherZ}

\title{
  Fisher's Z transformation 
}

\description{
  Fisher's Z transformation on an estimator in the range of 0~1. 
}
\usage{ FisherZ(x)}
\arguments{
  \item{x}{A numeric value, an estimator within the range of 0 to 1.}
}
\details{ Fisher's Z transformation transform an estimator x which is in
  the range of 0 and 1 to
  \eqn{z=\frac{1}{2}\log(\frac{(1+x)}{(1-x)})}{z=1/2*log((1+x)/(1-x))}. The    
 transformation is useful for skewed distribution of a variable between
 0 and 1. z after transformation will be more normally distributed. Fisher's Z
 transformation has been implemented on correlation coefficient, volume
 under ROC. }

\value{
  Return a numeric value---the transformed estimator z.
}

\references{

  Fisher, R.A. (1915) Frequency distribution of the values of the
  correlation coefficient in samples of an indefinitely large
  population.  \emph{Biometrika} \bold{10} (4) 507
  
  Xiong, C. and van Belle, G. and Miller, J.P. and Morris, J.C. (2006)
  Measuring and Estimating Diagnostic Accuracy When There 
  Are Three Ordinal Diagnostic Groups. \emph{Statistics In Medicine} \bold{25}
  7 1251--1273.

}
\author{
 Jingqin Luo
}

\note{
 Bug reports, malfunctioning, or suggestions for further improvements or
 contributions can be sent to Jingqin Luo <rosy@wubios.wustl.edu>.
}

\seealso{
 \code{\link{FisherZ.Var}}
}
\examples{

x <- 0.9
z <- FisherZ(x)
}
\keyword{manip}
\keyword{misc}

