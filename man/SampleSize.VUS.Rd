\name{SampleSize.VUS}
\alias{SampleSize.VUS}

\title{
  Sample size for volume under surface (VUS) 
}
\description{
Calculates the sample size required to estimate volume under surface
(VUS) within a given margin of error when the diagnostic marker follows
normal distributions in each diagnosis group.  
}
\usage{
SampleSize.VUS(mu.minus, mu0, mu.plus, s.minus, s0, s.plus, p = 0, q =
0, lam.minus = 1/3, lam0 = 1/3, lam.plus = 1/3, typeIerror = 0.05,
 margin = 0.05, subdivisions = 50000, ...)  
}

\arguments{
  \item{mu.minus}{
 A numeric value, the sample mean estimate for the mean parameter of normal
    distribution in \eqn{D^-}{D-}.  
}
  \item{mu0}{
 A numeric value, the sample mean estimate for the mean parameter of normal
  distribution in \eqn{D^0}{D0}.
}
  \item{mu.plus}{
 A numeric value, the sample mean estimate for the mean parameter of normal
  distribution in \eqn{D^+}{D+}.
}
  \item{s.minus}{
 A numeric value, the sample standard deviation (SD) estimate for the SD
  parameter of normal distribution in \eqn{D^-}{D-}.
}
  \item{s0}{
 A numeric value, the sample SD estimate for the SD parameter of normal
  distribution in \eqn{D^0}{D0}.
}
\item{s.plus}{
 A numeric value, the sample SD estimate for the SD parameter of normal
  distribution in \eqn{D^+}{D+}.
}
  \item{p}{
 A numeric value, the minimum required specificity, \eqn{0\le p < 1}{0<=p<1}, for calculation
 partial volume under ROC surface.  Default, p=0. 
 
}
  \item{q}{
 A numeric value, the minimum desired sensitivity,\eqn{0\le q<1}{0<=q<1}, for calculation
  partial volume under ROC surface. Default, q=0.  p=q=0 will give the
  complete VUS estimate, otherwise give the partial VUS estimate
  satisfying specificity no less than p and sensitivity no less than q.    
}
  \item{lam.minus}{
 A numeric value, the expected population proportion of the \eqn{D_-} group, used for
  sample size calculation. Default, lam.minus=1/3. The proportions of the three ordinal
  groups (lam.minus, lam0, lam.plus) should sum to 1.    
}
  \item{lam0}{
 A numeric value, the expected population proportion of the \eqn{D_0} group, used for
  sample size calculation. Default, lam0=1/3.  The proportions of the three ordinal
  groups (lam.minus, lam0, lam.plus) should sum to 1.  
}
  \item{lam.plus}{
 A numeric value, the expected population proportion of the \eqn{D_+} group, used for
  sample size calculation. Default, lam.plus=1/3.The proportions of the three ordinal
  groups (lam.minus, lam0, lam.plus) should sum to 1.    
}
  \item{typeIerror}{
 A numeric value, (1-typeIerror)*100\% confidence interval (CI) in sample size
  calculation. Default typeIerror=0.05, i.e., calculate 95\% CI.  
}
\item{margin}{
 A numeric value, the margin of error on the VUS estimates in sample size
  calculation. Default, margin=0.05. 
}
  \item{subdivisions}{
    The number of subintervals for integration using
    adaptive quadrature in the R function
    \code{\link{integrate}}. Default, subdivisions=50000.   
}
\item{\dots}{
  Other arguments that can be passed to the R function
    \code{\link{integrate}}, e.g., abs.tol, rel.tol, stop.on.error etc.     
}
}
\details{
  The (1-typeIerror)\% CI on VUS estimate under normality is
  \eqn{(VUS-Z_a*SE(VUS),VUS+Z_a*SE(VUS))}, then the sample size will be
  calculated such that : margin= \eqn{Z_a*SE(VUS)} or half of
  the CI's length, where \eqn{Z_a} is the normal quantile, \eqn{Z_a=1.95} given
  default typeIerror a=0.05.  
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
 \code{\link{Normal.VUS}}
}
\examples{

 mu.minus <- 6
 mu0 <- 8
 mu.plus <- 13.36


 s.minus <- 1
 s0 <- 1.2
 s.plus <- 1.4

 n <- SampleSize.VUS(mu.minus,mu0,mu.plus,s.minus,s0,s.plus,p=0,q=0,typeIerror=0.05) 
}
\keyword{design}
\keyword{univar}

