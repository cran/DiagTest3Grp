\name{Youden3Grp.OptimalCutoff.Normal.trueMeanSD}
\alias{Youden3Grp.OptimalCutoff.Normal.trueMeanSD}

\title{
  The Youden index and optimal cut point estimation under normal
  distribution for three ordinal diagnostic groups
}
\description{
 For a diagnostic test folling normal distributions, given the true or
 estimated normal distribution mean and standard deviation (SD) parameters for three
 ordinal diagnostic groups, ,calculate the Youden index and optimal cut
 point estimations for three ordinal diagnostic groups (usllay healthy
 D- , mildly diseased D0 and diseased D+).
}
\usage{ Youden3Grp.OptimalCutoff.Normal.trueMeanSD(mu.minus, mu0,
 mu.plus, s.minus, s0, s.plus) 
}

\arguments{
  \item{mu.minus}{
A numeric value. The mean parameter of normal distribution in D- group.
}
\item{mu0}{
A numeric value. The mean parameter of normal distribution in D0 group.
}
\item{mu.plus}{
 A numeric value. The mean parameter of normal distribution in D+ group.
}
\item{s.minus}{
A numeric value.  The SD parameter of normal distribution in D- group.
}
\item{s0}{
A numeric value. The SD parameter of normal distribution in D0 group.
}
\item{s.plus}{
A numeric value. The SD parameter of normal distribution in D+ group.
}
}
\details{
 See references.
}
\value{
  A data frame of numeric components,
  \item{t.minus}{The dervied optimal lower cut-point t.minus}
  \item{t.plus}{The dervied optimal upper cut-point t.minus}
  \item{youden}{The dervied optimal youden index}
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
 \code{\link{Youden3Grp}} 
}
\examples{

 mu.minus <- 6
 s.minus <- 1

 mu0 <- 8
 s0 <- 1.2

 mu.plus <- 11.05
 s.plus <- 1.4


Youden3Grp.OptimalCutoff.Normal.trueMeanSD(mu.minus,mu0,mu.plus,s.minus,s0,s.plus)

}
\keyword{htest}

