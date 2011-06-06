\name{NonParametric.VUS}
\alias{NonParametric.VUS}

\title{
 Nonparametric estimate of volumn under ROC surface (VUS)
}
\description{
   Nonparametric estimate of volumn under ROC surface (VUS)  
}
\usage{
 NonParametric.VUS(x, y, z,alpha=0.05,NBOOT=50)  
}

\arguments{
  \item{x}{
   A numeric vector, a diagnostic test's measurements in  the \eqn{D^-} (usually healthy subjects).
}
\item{y}{
 A numeric vector, a diagnostic test's measurements in  the \eqn{D^0} (usually mildly diseased subjects).
 }
 \item{z}{
   A numeric vector, a diagnostic test's measurements in  the \eqn{D^+}
   (usually severely diseased subjects).   
 }
 \item{alpha}{A numeric value, (1-alpha)*100\% percentile Confidence interval of the VUS estimate
   under normal assumption. Default, alpha=0.05.  }
 \item{NBOOT}{A integer.Total number of bootstrap samples for implementation.}

 }
\details{
The volume under ROC surface (VUS) indicates the probability of
correctly ranking a randomly selected  triplet \eqn{(U,V,W)} of
a diagnostic test's measurements, each from
\eqn{D^-},\eqn{D^0},\eqn{D^+}, i.e., \eqn{VUS=Pr\{U<V<W\}}. The
nonparametric estimator estimates VUS by empirical CDF method,
calculating the proportion of correct orderings among all possible
triplets, each from three diagnosis groups.    
}
\value{
  Return a numeric value as the nonparametric estimate of the VUS.
}
\references{
  Xiong, C. and van Belle, G. and Miller, J.P. and Morris, J.C. (2006)
  Measuring and Estimating Diagnostic Accuracy When There 
  Are Three Ordinal Diagnostic Groups. \emph{Statistics In Medicine} \bold{25}
  7 1251--1273.

  Ferri, C. and Hernandez-Orallo, J. and Salido, M.A. (2003) Volume
  under the ROC Surface for Multi-class Problems \emph{LECTURE NOTES IN
    COMPUTER SCIENCE} 108--120. 
}

\author{
  Jingqin Luo
}
\note{
 Bug reports, malfunctioning, or suggestions for further improvements or
 contributions can be sent to Jingqin Luo <rosy@wubios.wustl.edu>.   
}

\seealso{
 \code{\link{VUS}}  \code{\link{Normal.VUS}} \code{\link{NonParametric.VUS.var}}
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

 nonpar.vus <- NonParametric.VUS(x=x,y=y,z=z)
}

\keyword{nonparametric}
\keyword{htest}

