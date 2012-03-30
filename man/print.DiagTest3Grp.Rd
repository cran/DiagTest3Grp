\name{print.DiagTest3Grp}
\alias{print.DiagTest3Grp}
\title{Print an object of class DiagTest3Grp  
}
\description{
S3 class method print for a \code{\link{DiagTest3Grp}} object: Print
relevant information on a DiagTest3Grp object.
}
\usage{
\method{print}{DiagTest3Grp}(x,...)
}

\arguments{
  \item{x}{an object of class DiagTest3Grp}
  \item{...}{Other arguments for \code{\link{plot}} function}
}
\value{
  Print: the type of the three-group diagnostic test summary measure
  (VUS or Youden); the method used for the estimation; the data summary
  on group sample size, mean and SD; summary measure estimate and associated
  Confidence interval (CI); optimal cut-points; the three correct
  classification probability for each diagnosis group; and sample size
  (if using "Normal" method) to estimate a summary measure within
  specified precision. 
}

\author{
   Jingqin Luo
}

\seealso{
 \code{\link{VUS}} \code{\link{Youden3Grp}}
}
\keyword{print}
