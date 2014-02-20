\name{VUS}
\alias{VUS}

\title{
 A wrapper function for Volume under Surface (VUS) estimate, variance
 estimate under normal and nonparametric assumption and sample size calculation
}
\description{
   A wrapper function to calculate the Volume under Surface (VUS)
   estimate, its variance estimate and optimal cut-point, under normal
   and nonparametric assumption, to provide partial VUS estimate with a minimum
   requirement on the specificity and sensitivity under normality and to
   calculate the sample size under normality to achieve a certain
   estimation precision on VUS estimate. }
\usage{ VUS(x, y, z, method =c("Normal","NonPar"), p = 0, q = 0, alpha = 0.05,
          NBOOT = 100, subdivisions = 50000, lam.minus = 1/3, lam0 =1/3,
          lam.plus= 1/3, typeIerror = 0.05, margin =0.05,FisherZ=FALSE,
          optimalCut=TRUE,cut.seq=NULL,optimize=FALSE,...) 

}

\arguments{
  \item{x}{
    A numeric vector, a diagnostic test's measurements in  the D- (usually healthy subjects).
}
  \item{y}{
A numeric vector, a diagnostic test's measurements in  the D0 (usually mildly diseased subjects).
}
  \item{z}{
A numeric vector, a diagnostic test's measurements in  the D+ (usually severely diseased subjects).
}
  \item{method}{a character argument, method =\dQuote{Normal} or \dQuote{NonPar}, estimate the VUS under
    normality or nonparametrically.
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
A numeric value, (1-alpha)*100\% Confidence interval of the VUS estimate
  under normal assumption. Default, alpha=0.05.  
}
  \item{NBOOT}{ A numeric value, the total number of bootstrapping used to estimate the
    variance of nonparametric estimate of VUS. 
}
\item{subdivisions}{
A numeric value, the number of subintervals for integration using
  adaptive quadrature in the R function
  \code{\link{integrate}}. Default, subdivisions=50000.  
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
    calculation. Default, margin=0.05. The (1-typeIerror)\% CI on VUS
    estimate under normality is
    \eqn{(VUS-Z_a*SE(VUS),VUS+Z_a*SE(VUS))}, then margin=\eqn{Z_a*SE(VUS)} or half of
    the CI's length, where Z_a is the normal quantile, \eqn{Z_a=1.95} given
    default typeIerror a=0.05.
  }
  \item{FisherZ}{Reference to the argument in \code{\link{Normal.VUS}}.}
  \item{optimalCut}{
   A logic value of TRUE or FALSE. If TRUE, the
    function will return optimal cut-point from VUS analyses.
  }
  \item{cut.seq}{A sequence of numeric values from which the optimal cut-point
    will be selected from, by default=NULL, will use the unique values
    of the collection of x,y,z. }
  \item{optimize}{A logical value of TRUE or FALSE. If FALSE, take the
    empirical optimal cut point identified by empirical search within the given
    cut.seq as final reported optimal cut point.If TRUE, using the
    empirical optimal cut point as starting point in optimization algorithm for final
    optimal cut point.}
  \item{...}{
  Other arguments that can be passed to the R function
  \code{\link{integrate}}, e.g., abs.tol, rel.tol, stop.on.error etc. }
}
\details{
 For three ordinal group diagnostic test, there are two underlying
cut-point \eqn{t_-} and \eqn{t_+} with \eqn{t_-<t_+} based on which
patients are divided into the three ordinal groups. Patients
with a diagnostic test below \eqn{t_-} will be assigned to \eqn{D^-};those
with the test above \eqn{t_+} will be assigned to \eqn{D^+} and the remaining fall
into \eqn{D^0}. Following the specificifity and sensitity definition as
in diagnostic test for two groups, we call the probability of the first
two events as specificity \eqn{x=P_-\{T \le t_-\}=F_-(t_-)} and
sensitivity.\eqn{y=P_+\{T > t_+\}=1-F_+(t_+)=G_+(t_+)} where the
\eqn{P_i} and \eqn{F_+} denotes the probability density function and
cumulative density function of a diagnostic test in
\eqn{D^i}, \eqn{i=-,0,+} separately. Then, the probability that a
patient randomly selected from the \eqn{D^0} group has the test result
between the two cut-points can be expressed as,  
\eqn{z=P_0\{t_- \le T \le t_+\}=F_0(t_+)-F_0(t_-)=F_0(G_+^{-1}(y))-F_0(F_-^{-1}(x))} 

where the notation \eqn{H^{-1}(.)}{H^{-1}(.)} denotes the inverse function of \eqn{H}{H}
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
G_+(F_-^{-1}(x))\}}. The equation of partial VUS will be similar to the
above but the integration domain is \eqn{D_{pq}=\{p\le x \le 1,q \le y \le
  G_+(F_-^{-1}(x))\}}.  
The optimal cut-points from VUS analyse are defined as the one
}
\value{
  A object of \code{\link{DiagTest3Grp}} with a list of components. 

  \item{type}{A character value, type=dQuote{VUS} for VUS and type=dQuote{Youden} for the extended Youden index, indicating which summary measure is outputted.}
  \item{method}{A character value. For VUS, method can be
    \dQuote{Normal} or \dQuote{NonPar} (nonparametric); for Youden
    index, choices are \dQuote{Normal/TN/EMP/KS/KS-SJ}, indicating which
    method is used to estimate the summary measure.}
   \item{dat}{A list of 3 components. Three components have names \dQuote{x},\dQuote{y},\dQuote{z}, each recording the
     inputted marker measurements (after removing NAs) under \eqn{D^-},\eqn{D^0},\eqn{D^+} respectively.}
   \item{dat.summary}{A data frame with 3 rows (D-, D0,D+) and 3 columns (number of observations,mean, SD).}
   \item{estimate}{A numerical value. Point estimate for the summary
     measure, either VUS or Youden.}
   \item{variance}{A numeric value. Variance on the summary measure
     estimate. For normal method, output normal variance; for other methods output variance from bootstrapping.}    
   \item{CI}{A named numeric vector of length 2. confidence interval on
     the summary measure estimate, with name like 2.5\%, 97.5\% if
     significance level is set to be 5\%. For both VUS and the Youden
     index, when normal method is in use, the CI is normal CI while
     bootstrap method was used under other methods.}
   \item{cut.point}{A named numeric vector of length 2. optimal
     cut-points with name \dQuote{t.minus} for lower optimal cut point
     and name \dQuote{t.plus} for upper optimal point.}
   \item{classify.prob}{ A named numeric vector of 3 values. Estimates on the
     three group correct classification probabilities. specificity on \eqn{D^-}:
     \eqn{Sp==Pr(x\le t_-|D^-)}; sensitivity on \eqn{D^+}:
     \eqn{Se=Pr(z\ge t_+|D^+)}; correct classification probability on
     \eqn{D^0}: \eqn{Sm=SPr(t_-<y<t_+|D^0)}. For VUS, it's empirical
     estimation. For Youden index, depending on method adopted for the
     Youden index estimate, the three probabilities will be estimated
     using specified method.} 
   \item{sampleSize}{A numeric value The sample size to estimate the summary measure within given margin of error and type-I error rate.See \code{\link{SampleSize.VUS}} and \code{\link{SampleSize.Youden3Grp}}.}
   \item{alpha}{A numeric value. The significance level for the CI
     computation, e.g, default=5\%.}
   \item{typeIerror}{A numeric value for type-I error rate, e.g.,default=5\%.}
   \item{margin}{A numeric value. The margin of errors (precision) to
     estimate the summary measure s.t. the half the length of the
     resulting CI is equal to the given margin.}
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

  Ferri, C. and Hernandez-Orallo, J. and Salido, M.A. (2003) Volume
  under the ROC Surface for Multi-class Problems \emph{LECTURE NOTES IN
    COMPUTER SCIENCE} 108--120. 

}
\author{
 Jingqin Luo
}
\note{
  Bug reports, malfunctioning, or suggestions for further improvements or
  contributions can be sent to Jingqin Luo <rosy@wubios.wustl.edu>.  
}

\section{Warning }{The bootstrapping to obtain the variance on the
  nonparametric VUS estimate may take a while. }

\seealso{
 \code{\link{Normal.VUS}} \code{\link{NonParametric.VUS}}  \code{\link{NonParametric.VUS.var}}
}
\examples{

 data(AL)
 group <- AL$group
 table(group)

 ##take the negated kfront marker measurements
 kfront <- -AL$kfront

 x <- kfront[group=="D-"]
 y <- kfront[group=="D0"]
 z <- kfront[group=="D+"]

 ##normal estimate
 normal.res <- VUS(x,y,z,method="Normal",p=0,q=0,alpha=0.05)
 normal.res

 ##nonparametric estimate
\dontrun{
 nonpar.res <- VUS(x,y,z,method="NonPar",p=0,q=0,alpha=0.05,NBOOT=100)
 nonpar.res
}

 ## S3 method for class 'DiagTest3Grp':
 print(normal.res)

 ## S3 method for class 'DiagTest3Grp':
 plot(normal.res)

}
\keyword{htest}
\keyword{nonparametric}
\keyword{design}
