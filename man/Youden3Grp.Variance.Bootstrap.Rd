\name{Youden3Grp.Variance.Bootstrap}
\alias{Youden3Grp.Variance.Bootstrap}

\title{
 Obtain bootstrap variance and confidence itnerval on the extended
 Youden index and optimal cut-off points  }

\description{
Given the diagnostic test measurements x, y, z in the three ordinal groups
\eqn{D^-}, \eqn{D^0}, \eqn{D^+} separately, the function repeatedly draws a bootstrap sample
each of x, y, z to estimate the extended Youden index for three ordinal
groups and the associated optimal cut point and obtain the final
bootstrap average estimate on the Youden index and optimal cut-off point
and the confidence interval.
}
\usage{
Youden3Grp.Variance.Bootstrap(x, y, z, method = "Normal", seed.seq =NULL,
randomStart.N = 5, NBOOT=10,alpha=0.05,t.minus.start = NULL,
t.plus.start = NULL, ...)  
}

\arguments{
  \item{x}{A numeric vector, a diagnostic test's measurements for
    subjects in D- (usually a healthy group).    
}
\item{y}{A numeric vector, a diagnostic test's measurements for subjects in D0
    (usually a mildly diseased group).  
}
\item{z}{A numeric vector, a diagnostic test's measurements for subjects in D+
    (usually a severly diseased group).  
}
  \item{method}{
A character argument. Specify a method to be used for estimating the extended Youden index
    Choices include
    \sQuote{Normal}---Estimate the extended Youden index under the
      assumption of  normal distributions of a diagnostic test in
      the three groups.
    
    \sQuote{TN}---Transformed normal. Implement Box-cox transformation to
      approximate normality and then estimate the extended Youden index
      under normality.
    
    \sQuote{EMP}---Estimate the extended Youden index by using empirical
      cumulative density function.
 
    \sQuote{KS}---Estimate the extended Youden index by using Kernel
      density estimation with a normal reference rule for bandwidth
      selection.
    

    \sQuote{KS-SJ}---Estimate the extended Youden index by using Kernel
      density estimation with the Sheather-Jones Plug-in method for
      bandwidth selection.
  }  
 
  \item{seed.seq}{A numeric vector. Users can specify a sequence of random seeds for
    bootstrap sampling of x. Boostrapping of y and z will use seed0+1
    and seed0+2. Default will be the sequence of 1:10.  }
  
  \item{randomStart.N}{A numeric value. An argument need to be specified when the method
    \dQuote{EMP}, \dQuote{KS} or \dQuote{KS-SJ} is used to estimate the
    optimal cut-point and the extended Youden index J. Default,
    randomStart.N=1.  See \code{\link{Youden3Grp}}.  
  }
  \item{NBOOT}{A numeric value. Total number of bootstrapping, default=10.}
  \item{t.minus.start}{ The starting points of the lower optimal
    cut-point (t-) which separate the D- and D0 group. Default
    t.minus.start=NULL will randomly generate starting value. See
    \code{\link{Youden3Grp}}.   
  }

  \item{t.plus.start}{
A numeric value. The starting points of the upper optimal cut-point (t+) which separate the
    D0 and D+ group. Default t.plus.start=NULL will randomly generate
    starting value. See \code{\link{Youden3Grp}}.
  }
  \item{alpha}{A numeric value. the significance level, will provide the basic quantile
    confidence interval (alpha/2*100\%,1-alpha/2*100\%).}
  \item{\dots}{
    Other arguments to be passed to the R function \code{\link{Youden3Grp}}.
  }
}
\details{
This function is carried to get bootstrap estimates of the extended
Youden index and associated cut-points to provide confidence interval. 
See details in \code{\link{Youden3Grp}}
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
}

\references{
Luo, J and Xiong, C. 2010 Youden Index and Associated Optimal Cut-point
for Three Ordinal Groups. \emph{Communications In Statistics-Simulation
  and Computation}.
}

\author{  
  Jingqin Luo
}

\note{
 Bug reports, malfunctioning, or suggestions for further improvements or
 contributions can be sent to Jingqin Luo <rosy@wubios.wustl.edu>.      
}

\seealso{
  \code{\link{Youden3Grp}}  \code{\link{BW.ref}}
}

\examples{
x <- rnorm(50,6,1)

y <- rnorm(60,8,1.2)

z <- rnorm(40,10,1.4)

temp.res <- Youden3Grp.Variance.Bootstrap(x=x, y=y, z=z,
method="Normal",seed.seq=1:10,randomStart.N=3,NBOOT=10)  
}

\keyword{univar}
\keyword{htest}
