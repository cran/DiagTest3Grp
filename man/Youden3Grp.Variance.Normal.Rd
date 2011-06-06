\name{Youden3Grp.Variance.Normal}
\alias{Youden3Grp.Variance.Normal}

\title{
 Variance estimations for the extended Youden index and associated lower
 and upper optimal cut-point under normality assumption
}
\description{
 Calculate the variance estimation for the extended Youden index and the
 variance estimation for the associated lower
 and upper optimal cut-point assuming that a diagnostic test follows
 normal distributions in the three ordinal groups
 (\eqn{D^-},\eqn{D^0},\eqn{D^+}).  
}
\usage{
  Youden3Grp.Variance.Normal(x, y, z, alpha = 0.05)
}

\arguments{
  \item{x}{
A numeric vector. A diagnostic test's measurements in  the \eqn{D_-}{D-} (usually healthy
    subjects).
  }

  \item{y}{
A numeric vector. A diagnostic test's measurements in  the \eqn{D^0}{D0} (usually mildly
  diseased subjects). 
}
\item{z}{
A numeric vector. A diagnostic test's measurements in  the \eqn{D^+}{D+} (usually severely
   diseased subjects).  
}
\item{alpha}{
A numeric value. Significance level so that the function calculates the (1-alpha)*100\% confidence interval (CI) on the estimates of
  the extended Youden index and optimal cut-points under normal
  assumption. Default, alpha=0.05.   
 }
}
\details{
  See details in \code{\link{Youden3Grp}}.
}
\value{ 
  Return a list, including the following components
  \item{var.youden}{The normal-method based variance on the optimal Youden index.}
  
  \item{var.t.minus}{A numeric value.The variance on the lower optimal cut-point t.minus.}

  \item{var.t.plus}{A numeric value.The variance on the upper optimal cut-point t.plus.}

  \item{var.youden.z}{A numeric value. The variance on the Fisher's Z transformed optimal Youden index.}

  \item{youden.CI}{A named numeric of length 2. CI for the estimate of youden with the lower and the upper
    CI limit (under the names alpha/2*100\% and (1-alpha/2)*100\%).}

  \item{t.minus.CI}{A named numeric of length 2. CI for the estimate of t.minus  (the lower
    optimal cut-point) with the lower and the upper
    CI limit (under the names alpha/2*100\% and (1-alpha/2)*100\%) CI}
  \item{t.plus.CI}{A named numeric of length 2. CI for the estimate of t.plus (the upper
    optimal cut-point) with the lower and the upper
    CI limit (under the names alpha/2*100\% and (1-alpha/2)*100\%).}
  \item{youden.z.CI}{A named numeric of length 2. CI for the estimate of
    Fisher-Z transformed youden with the lower and the upper
    CI limit (under the names alpha/2*100\% and (1-alpha/2)*100\%).}
   \item{partialDeriv}{A numeric data frame with one row and multiple
     columns, containing estimated SD parameters in each diagnosis group and the partial derivatives of Youden estimate
     w.r.t the relevant mean and SD parameters which are outputted for performance
     of statistical tests on markers under normal method or NA under
     nonparametric method.}
}
\references{
  Luo, J and Xiong, C. 2010 Youden Index and Associated Optimal
  Cut-point for Three Ordinal Groups. \emph{Communications In
    Statistics-Simulation and Computation}.
}

\author{
 Jingqin Luo
}

\note{
  Bug reports, malfunctioning, or suggestions for further improvements or
 contributions can be sent to Jingqin Luo <rosy@wubios.wustl.edu>.    
}

\seealso{
 \code{\link{Youden3Grp}} \code{\link{Youden3Grp.Variance.Bootstrap}}
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


 temp.res <- Youden3Grp.Variance.Normal(x=x, y=y, z=z, alpha=0.05)

 ###variance of the extended Youden index and optimal cut-points
 temp.res[c("var.youden","var.t.minus","var.t.plus")]

  }

\keyword{ univar }
\keyword{ htest}
