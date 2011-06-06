\name{FisherZ.Var}
\alias{FisherZ.Var}

\title{
 Calculate the variance of an estimator after Fisher's Z transformation
}
\description{
 Calculate the variance of an estimator after Fisher's Z transformation

}
\usage{
  FisherZ.Var(x, var.x)
}

\arguments{
  \item{x}{
  A numeric value, an estimator in the range of 0~1.
}
  \item{var.x}{
     A numeric value, the variance of the estimator x.
}
}
\details{
 The variance of an estimator x after Fisher's Z transformation is  \eqn{\frac{var(x)}{(1-x^2)^2}{var(x)/(1-x^2)^2}}.
}
\value{
 Return a numerical value---the variance of an estimator after Fisher's Z transformation.
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

\seealso{
 \code{\link{FisherZ}}
}
\examples{

x <- 0.8

var.x <- 0.06

new.var <- FisherZ.Var(x,var.x)
}
\keyword{univar}
\keyword{manip}
\keyword{misc}
