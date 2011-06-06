\name{bootSample}
\alias{bootSample}
\title{Draw a bootstrap sample from a given vector of data values  
}
\description{
Draw one bootstrap sample of a given vector of sample values    
}
\usage{
  bootSample(x, seed0)
}

\arguments{
  \item{x}{a numeric vector, sample data.}
  \item{seed0}{a random seed for bootstrap sampling.}
}
\value{
  Returns a numeric vector which is a bootstrap sample drawn randomly
  with replacement from and with the same vector length as the sample data  x .   
}

\author{
   Jingqin Luo
}

\examples{

##generate random data
x <- rnorm(100)

##draw a bootstrap sample
boot.x <- bootSample(100,seed0=12345)
}
\keyword{datagen}
