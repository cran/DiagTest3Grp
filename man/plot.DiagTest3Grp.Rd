\name{plot.DiagTest3Grp}
\alias{plot.DiagTest3Grp}
\title{Scatter plot and boxplot on an object of class DiagTest3Grp}
\description{
  S3 class method plot for a \code{\link{DiagTest3Grp}} object:
  Scatter-plot and box-plot}

\usage{
 plot.DiagTest3Grp(x,...)
}

\arguments{
  \item{x}{an object of class DiagTest3Grp}
  \item{...}{Other arguments for \code{\link{plot}} function}
}
\value{
Scatter-plot and Box-plot of the data points, with green, blue, red
color corresponding to group \eqn{D^-}, \eqn{D^0}, \eqn{D^+} and summary
measure estimate as well as optimal cut points are indicated.
}

\author{
   Jingqin Luo
}

\seealso{
 \code{\link{VUS}} \code{\link{Youden3Grp}}
}

\keyword{hplot}
