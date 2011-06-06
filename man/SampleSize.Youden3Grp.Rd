\name{SampleSize.Youden3Grp}
\alias{SampleSize.Youden3Grp}

\title{
  Sample size for the extended Youden index for three ordinal groups
}
\description{
  Calculates the sample size required to estimate  the extended Youden
  index within a given margin of error when the diagnostic marker follows
  normal distributions in each diagnosis group. 
}
\usage{SampleSize.Youden3Grp(mu.minus, mu0, mu.plus, s.minus, s0, s.plus,
t.minus, t.plus, lam.minus, lam0, lam.plus, typeIerror = 0.05,
 margin=0.05)  
}

\arguments{
  \item{mu.minus}{ A numeric value, the sample mean estimate for the mean parameter of normal
    distribution in \eqn{D^-}{D-}. }
  \item{mu0}{A numeric value, the sample mean estimate for the mean parameter of normal
  distribution in \eqn{D^0}{D0}.  }
  \item{mu.plus}{A numeric value, the sample mean estimate for the mean parameter of normal
  distribution in \eqn{D^+}{D+}.  }
  \item{s.minus}{A numeric value, the sample standard deviation (SD) estimate for the SD
    parameter of normal distribution in \eqn{D^-}{D-}. }
  \item{s0}{A numeric value, the sample SD estimate for the SD parameter of normal
  distribution in \eqn{D^0}{D0}.}
  \item{s.plus}{A numeric value, the sample SD estimate for the SD parameter of normal
  distribution in \eqn{D^+}{D+}.}
  \item{t.minus}{A numeric value, the optimal lower cut point to attain the optimal
    extended Youden index under normal assumptions.}
 \item{t.plus}{A numeric value, the optimal upper cut point to attain the optimal
   extended Youden index under normal assumptions.   }
 \item{lam.minus}{A numeric value, the expected population proportion of the \eqn{D_-}
   group, used for sample size calculation. Default, lam.minus=1/3. The
   proportions of the three ordinal groups (lam.minus, lam0, lam.plus)
   should sum to 1.  } 
 \item{lam0}{A numeric value, the expected population proportion of the \eqn{D_0} group,
   used for sample size calculation. Default, lam0=1/3.  The proportions of the
  three ordinal groups (lam.minus, lam0, lam.plus) should sum to 1. }
  \item{lam.plus}{A numeric value, the expected population proportion of the \eqn{D_+}
    group, used for sample size calculation. Default, lam.plus=1/3.The
    proportions of the three ordinal groups (lam.minus, lam0, lam.plus)
    should sum to 1. } 
  \item{typeIerror}{A numeric value, (1-typeIerror)*100\% confidence interval
    (CI) based on which the length of CI for sample size  calculation
    will be obtained. Default typeIerror=0.05, i.e., calculate  95\% CI. }
  \item{margin}{A numeric value, the margin of error on the extended Youden index
    estimates in sample size calculation. Default, margin=0.05. }
}

\details{
  The (1-typeIerror)\% CI on the Youden index (J) estimate under normality is
  \eqn{(J-Z_a*SE(J),J+Z_a*SE(J))}, then margin= \eqn{Z_a*SE(J)} or half of
  the CI's length, where \eqn{Z_a} is the normal quantile, \eqn{Z_a=1.95} given
  default typeIerror \eqn{a=0.05}.  
}

\value{
   Returns a numeric value indicating required sample size.
}

\references{
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
 \code{\link{Youden3Grp}}
}
\examples{

 mu.minus <- 6
 mu0 <- 8
 mu.plus <- 13.36


 s.minus <- 1
 s0 <- 1.2
 s.plus <- 1.4

 t.minus <- 7.02
 t.plus <- 10.52

 n <- SampleSize.Youden3Grp(mu.minus, mu0, mu.plus, s.minus, s0, s.plus,
t.minus, t.plus)

}
\keyword{design}
\keyword{univar}

