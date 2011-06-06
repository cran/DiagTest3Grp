\name{var.mu}
\alias{var.mu}
  
\title{
  Variance estimate of the sample estimate for a normal mean parameter  
}
\description{
  Calculates the variance estimate of a sample estimate for a normal mean
  parameter.   
}
\usage{
var.mu(s, n)
}

\arguments{
  \item{s}{
    A numeric value, the sample estimate for normal standard deviation.
}
  \item{n}{
A numeric value, the sample size of data.
}
}
\details{
 For a sample dataset of size n drawn from a normal distributin
 \eqn{N(\mu,\sigma^2)}{N(\mu,\sigma^2}, denote the sample estimate on 
 the mean parameter and standard deviation parameter as \eqn{u}{u} and s,
 Then the variance of the sample mean is \eqn{Var(u)=\frac{s^2}{n}}{Var(u)=s^2/n}. 
}
\value{
  A numeric value, the variance of a sample estimate for a normal mean parameter.
}
\references{
  Patel, J.K. and Read, C.B. (1996) Handbook of the Normal
  Distribution. CRC Press.

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
\code{\link{var.sigma}}
}
\examples{
s <- 1.5
n <- 100
var.mu(s,n)
}
\keyword{htest}
\keyword{univar}
