\name{var.sigma}
\alias{var.sigma}

\title{ Variance of a sample estimate for a normal standard deviation
  (SD) parameter  } 
\description{ Calculates the variance of a sample estimate for a normal SD
  parameter.     } 
\usage{
 var.sigma(s,n)
 }
\arguments{
 \item{s}{A numeric value, the sample standard deviation estimate.}
 \item{n}{A numeric value, the sample size of data.}
}

\details{
 For a sample dataset of size n drawn from a normal distributin
 \eqn{N(\mu,\sigma^2)}{N(\mu,\sigma^2)}, denote the sample estimate on 
 standard deviation (SD) parameter as s. Then the variance of the sample
 SD estimate s is \eqn{Var(s)=\frac{s^2}{2(n-1)}}{Var(s)=0.5*s^2/(n-1)}.  }
\value{
  A numeric value, the variance estimate of the sample estimate for a normal SD parameter.
}
\references{
  Patel, J.K. and Read, C.B. (1996) Handbook of the Normal
  Distribution. CRC Press.
  
  Luo, J and Xiong, C. (2010) Youden Index and Associated Optimal Cut-point
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
\code{\link{var.sigma}}
}
\examples{
 s <- 1.5
 n <- 100
 var.sigma(s,n)
}
\keyword{ univar}
\keyword{ htest }
