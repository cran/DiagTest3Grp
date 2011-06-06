\name{Normal.VUS}
\alias{Normal.VUS}
  
\title{
 Volume under ROC surface (VUS) and partial VUS estimate and associated
 variance under normality assumptions
}
\description{
 Estimate volume under ROC surface (VUS), partial VUS and associated
 variance assuming a diagnostic test follows normal distribution in all
 three ordinal diagnosis groups \eqn{D_-}, 
 \eqn{D_0}, \eqn{D_+}.  
}
\usage{
 Normal.VUS(x, y, z, p = 0, q = 0, alpha = 0.05, subdivisions = 50000,
 lam.minus = 1/3, lam0 = 1/3, lam.plus = 1/3, typeIerror = 0.05,
 margin=0.05, ...)      
}

\arguments{
  \item{x}{
 A numeric vector, a diagnostic test's measurements in  the \eqn{D^-}{D-} (usually healthy
  subjects).    
}
\item{y}{
A numeric vector, a diagnostic test's measurements in  the \eqn{D^0}{D0} (usually mildly
  diseased subjects).  
}
  \item{z}{
A numeric vector, a diagnostic test's measurements in  the \eqn{D^+}{D+} (usually severely
   diseased subjects).  
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
\item{alpha}{
A numeric value, the (1-alpha)*100\% Confidence interval of the VUS estimate
  under normal assumption. Default, alpha=0.05.  
}
  \item{subdivisions}{A numeric value, the number of subintervals for integration using
    adaptive quadrature in the R function
    \code{\link{integrate}}. Default, subdivisions=50000. 
}
\item{lam.minus}{
A numeric value, the expected population proportion of the \eqn{D^-} group, used for
  sample size calculation. Default, lam.minus=1/3. The proportions of the three ordinal
  groups (lam.minus,lam0,lam.plus) should sum to 1.  
}
\item{lam0}{
A numeric value, the expected population proportion of the \eqn{D^0} group, used for
  sample size calculation. Default, lam0=1/3.
  The proportions of the three ordinal groups (lam.minus,lam0,lam.plus)
  should sum to 1.    
}
\item{lam.plus}{
A numeric value, the expected population proportion of the \eqn{D_+} group, used for
  sample size calculation. Default, lam.plus=1/3.The proportions of the
  three ordinal groups (lam.minus,lam0,lam.plus) should sum to 1.      
}
  \item{typeIerror}{
A numeric value, (1-typeIerror)*100\% confidence interval (CI) in sample size
    calculation. Default typeIerror=0.05, i.e., calculate 95\% CI.  
}
\item{margin}{A numeric value, the margin of error on the VUS estimates in sample size
  calculation. Default, margin=0.05. The (1-typeIerror)\% CI on VUS
  estimate under normality is
  \eqn{(VUS-Z_a*SE(VUS),VUS+Z_a*SE(VUS))}, then margin= Z_a*SE(VUS) or half of
  the CI's length, where \eqn{Z_a} is the normal quantile, Z_a=1.95 given
  default typeIerror=0.05.
}
\item{\dots}{Other arguments that can be passed to the R function
  \code{\link{integrate}}, e.g., abs.tol,rel.tol,stop.on.error etc.     
}
}
\details{
For three ordinal group diagnostic test, there are two underlying
cut-point \eqn{t_-} and \eqn{t_+} with \eqn{t_-<t_+} based on which
patients are divided into the three ordinal groups. Patients
with a diagnostic test below \eqn{t_-} will be assigned to \eqn{D^-}; those
with the test above \eqn{t_+} will be assigned to \eqn{D^+} and the remaining fall
into \eqn{D^0}. Following the specificifity and sensitity definition as
in diagnostic test for two groups, we call the probability of the first
two events as specificity \eqn{x=Pr_-\{T \le t_-\}=F_-(t_-)} and
sensitivity.\eqn{y=Pr_+\{T > t_+\}=1-F_+(t_+)=G_+(t_+)} where the
\eqn{Pr_i} and \eqn{F_+} denotes the probability density function (PDF) and
cumulative density function (CDF) of a diagnostic test in
\eqn{D^i}, \eqn{i=-,0,+} separately. Then, the probability that a
patient randomly selected from the \eqn{D^0} group has the test result
between the two cut-points can be expressed as,  
\eqn{z=P_0\{t_- \le T \le t_+\}=F_0(t_+)-F_0(t_-)=F_0(G_+^{-1}(y))-F_0(F_-^{-1}(x))} 

where the notation \eqn{H^{-1}(.)} denotes the inverse function of H.
z is a function of the specificity and sensitivity, i.e.,
\eqn{z=z(x,y)}, which constitutes a ROC surface in the three-dimensinal
space (x,y,z). The volume under the ROC surface (VUS) defined by z
can be written as,
\deqn{V_{00}=\int\int_{D_{00}} \{F_0(G_+^{-1}(y))-F_0(F_-^{-1}(x))\}d_x
  d_y=\int_{-\infty}^{+\infty}
  F_-(s)G_+(s)f_0(s)d_s}{V_{00}=\int\int_{D_{00}}
  \{F_0(G_+^{-1}(y))-F_0(F_-^{-1}(x))\}d_x d_y=\int_{-\infty}^{+\infty}
  F_-(s)G_+(s)f_0(s)d_s} 
The integration domain is \eqn{D_{00}=\{0\le x \le 1,0\le y \le
G_+(F_-^{-1}(x))\}}. Partial VUS requires a marker has specificity at
least p and sensitivity at least q. The equation of partial VUS will be similar to the
above but the integration domain is \eqn{D_{pq}=\{p\le x \le 1,q \le y \le
  G_+(F_-^{-1}(x))\}}.  
}
\value{ A list of components, including
  \item{dat}{A list of components \dQuote{x},\dQuote{y}, \dQuote{z},
    recording the original marker measurement inputted into the function.}
  \item{dat.summary}{A numeric data frame with three rows (\eqn{D^-},\eqn{D^0},
    \eqn{D^+}) and three columns (\dQuote{n},
    \dQuote{mu},\dQuote{sd}),summariz number of observations, mean and
    SD estimated for each diagnosis group from samples.}
  \item{estimate}{A numeric value.  A numeric value. If p=q=0, VUS gives the complete VUS estimate, otherwise give the partial VUS estimate
    satisfying specificity no less than p and sensitivity no less than
    q. }
  \item{variance}{A numeric value. The variance of VUS estimate under normal
    assumption.  }
  \item{CI}{A named numeric vector of length 2. The (1-alpha)*100\%  CI
    on VUS or partial VUS.}
  \item{sampleSize}{The sample size on investigating a diagnostic test
     based on the pilot measurments x, y, z. The sample size calculated such
     that based on given population proportion of the three diagnosis
     groups, the margin of error on the (1-typeIerror)*100\% confidence interval of the VUS 
     estimate under normality is within given margin of error (by
     default 0.05).}
   \item{partialDeriv}{A numeric data frame with one row and multiple
     columns, containing relevant
     parameters (a,b,c,d) and the partial derivatives of VUS estimate
     w.r.t the relevant parameters which are outputted for performance
     of statistical tests on markers under normal method or NA under
     nonparametric method.}
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
  \code{\link{VUS}}  \code{\link{NonParametric.VUS}}  \code{\link{NonParametric.VUS.var}}
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

 ##Complete VUS
 vus <- Normal.VUS(x=x, y=y, z=z)

 ##complete VUS estimate
 vus$VUS

 ##variance of the complete VUS estimate
 vus$variance

 ##sample size
 vus$sampleSize


 ##partial VUS with p>=0.4 and q>=0.5
 partial.vus <- Normal.VUS(x=x, y=y, z=z,p=0.4, q=0.5)

}
\keyword{htest}
\keyword{univar}
