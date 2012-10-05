\name{DiagTest3Grp-package}
\alias{DiagTest3Grp-package}
\alias{DiagTest3Grp}
\docType{package}
\title{
  Diagnostic test summary measures for three ordinal groups 
}
\description{
Assuming diagnostic marker increase monotonically and stochatically with
disease severity, the R package provides utilities to estimate two
diagnostic test summary measures for three ordinal groups ----volume
under ROC surface (VUS) and the extended Youden index. Variance,
confidence interval and optimal cut-points both under the normal
assumption and also the non-parametric method(s) will be provided for
the summary measures. Statistical tests are implemented to compare multiple diagnostic tests and
two diagnostic test. Sample size is calculated to estimate the summary
measure for a diagnostic test within user-specified margin of error for future study
planning.  
}

\details{
\tabular{ll}{
Package: \tab DiagTest3Grp\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-10-21\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr}


\code{\link{VUS}}

\code{\link{Youden3Grp}}

\code{\link{Youden3Grp.Variance.Bootstrap}}

\code{\link{SampleSize.VUS}}

\code{\link{SampleSize.Youden3Grp}}
}
\author{
 Jingqin Luo, Chengjie Xiong
 Maintainer: Jingqin Luo <rosy@wubios.wustl.edu>
}
\references{

  Xiong, C. and van Belle, G. and Miller, J.P. and Morris, J.C. (2006)
  Measuring and Estimating Diagnostic Accuracy When There 
  Are Three Ordinal Diagnostic Groups. \emph{Statistics In Medicine} \bold{25}
  7 1251--1273.

  Luo, J and Xiong, C. (2012) Youden Index and Associated Optimal Cut-point
  for Three Ordinal Groups. \emph{Communications In Statistics-Simulation
  and Computation}(in press).

  Jingqin Luo, Chengjie Xiong (2012). DiagTest3Grp: An R Package for
  Analyzing Diagnostic Tests with Three Ordinal Groups.\emph{Journal of
    Statistical Software}, \bold{51}(3), 1-24.URL
  http://www.jstatsoft.org/v51/i03/. 
}

\keyword{package}

\seealso{
  \code{\link{Youden3Grp}}
  \code{\link{VUS}}
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

 ## VUS under normality assumption
 normal.res <- VUS(x,y,z,method="Normal",p=0,q=0,alpha=0.05)
 ###S3 method
 print(normal.res)
 ###S3 method
 plot(normal.res)

 ##Youden index using normal method
 youden.res <- Youden3Grp(x=x,y=y,z=z,method="Normal")
 print(youden.res)
 plot(youden.res)

}
