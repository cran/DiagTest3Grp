\name{Youden3Grp}
\alias{Youden3Grp}
 
\title{
  The wrapper function for Youden index analysis
}
\description{
The wrapper function provide not only point estimates but also
variances for the extended Youden index for three ordinal
groups, the associated optimal cut-points under parametric and
non-parametric methods and sample size under normal.  
}
\usage{ Youden3Grp(x, y, z, method = c("Normal","TN","EMP","KS","KS-SJ"), randomStart.N = 1,
optim.method = NULL, t.minus.start = NULL, t.plus.start = NULL,
lam.minus = 1/3, lam0 = 1/3, lam.plus = 1/3, typeIerror = 0.05, margin
 =0.05,NBOOT=10,seed.seq=NULL,alpha=0.05, ...)  
}

\arguments{
  \item{x}{
    A numeric vector, a diagnostic test's measurements for
    subjects in D- (usually a healthy group).
}
  \item{y}{
    A numeric vector, a diagnostic test's measurements for subjects in D0
    (usually a mildly diseased group). 
}
\item{z}{
  A numeric vector, a diagnostic test's measurements for subjects in D+
    (usually a severly diseased group). 
}
  \item{method}{
    A character.
    Specify a method to be used for estimating the extended Youden index
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
    \item{randomStart.N}{ An argument need to be specified when the
      method \dQuote{EMP}, \dQuote{KS} or \dQuote{KS-SJ} is used to
      estimate the optimal cut-point and the extended Youden index
      J. Default,   randomStart.N=1.   }
  \item{optim.method}{A character argument. The optimization method specified to derive the
    optimal cut-points 
    when nonparametric methods are used to estimate the extended Youden
    index. Default, optim.method=L-BFGS-B.see \code{\link{optim}}. }
  \item{t.minus.start}{A numeric value. The starting points of the lower optimal
    cut-point (t-) which separate the D- and D0 group. Default
    t.minus.start=NULL will randomly generate  starting value. }
  \item{t.plus.start}{A numeric value. The starting points of the upper optimal
    cut-point (t+) which separate the 
 D0 and D+ group. Default t.plus.start=NULL will randomly generate
 starting value. }
\item{lam.minus}{A numeric value. The expected population proportion of the D- group,
  used for sample size calculation. Default, lam.minus=1/3. The
  proportions of the three ordinal groups (lam.minus, lam0, lam.plus)
  should sum to 1. }
  \item{lam0}{A numeric value. The expected population proportion of the D0 group, used for
  sample size calculation. Default, lam0=1/3. The proportions of the
  three ordinal  groups (lam.minus, lam0, lam.plus) should sum to 1. }
\item{lam.plus}{A numeric value. The expected population proportion of the D+ group, used for
  sample size calculation. Default, lam.plus=1/3.The proportions of the
  three ordinal groups (lam.minus, lam0, lam.plus) should sum to 1.}
  \item{typeIerror}{A numeric value. (1-typeIerror)*100\% confidence interval
    (CI) in sample size calculation. Default typeIerror=0.05, i.e.,
    calculate 95\% CI.  }  
  \item{margin}{A numeric value. The margin of error on the Youden index estimate in
    sample size calculation. Default, margin=0.05.  }
  \item{NBOOT}{A numeric value. Total number of bootstrap samples to
    draw for bootstrap variance and CI.}
  \item{seed.seq}{a numeric vector of length NBOOT.The sequence of random seeds for bootstrapping
    procedure to obtain bootstrapping variance and percentile confidence
    interval if methods other than \dQuote{Normal} is adopted}
  \item{alpha}{A numeric value. Significance level for CI
    calculation. Default, alpha=0.05}
  \item{...}{Other options for the R function \code{\link{optim}}.}
}

\details{
  Let the marker denoted by \eqn{T}{T} and the density function under
      \eqn{D^i}{Di} as \eqn{f_i}{fi}, \eqn{i=-,0,+}{i=-,0.+}. The Youden 
   index for three ordinal groups is defined as, 
  \eqn{J(t_-,t_+)=\frac{1}{2}[f_-(T\le t_-)+f_0(t_-\le T \le
    t_+)+f_+(T\ge t_+)-1)}{J(t-,t+)=1/2*(f_(T<= t-)+f0(t-<= T<= t+)+f+(T
    <= t+)-1]}. 
The optimal index achieved at an optimal pair will be adopted as the
final summary measure.  Note that the package assumes that the
 marker values increase monotonically and stochastically with disease
 severity (the marker means in the three group \eqn{D^-}{D-}, \eqn{D^0}{D0} and \eqn{D^+}{D+} are
 monotonically increasing). If a monotonically decreasing trend, the marker should be negated for anlaysis. 

 For the extended Youden index and optimal cut-point estimates, see
 reference paper for details.
 
 For sample size calculation,  the (1-typeIerror)\% CI on the Youden estimate under normality is
  \eqn{(J-Z_a*SE(J),J+Z_a*SE(J))}, then the sample size will be calculated
  such that margin=\eqn{Z_a*SE(J)} or half of
  the CI's length, where Z_a is the normal quantile, \eqn{Z_a=1.95} given
  default typeIerror \eqn{a=0.05}.
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
     resulting CI is equal to the given margin. 
     Note: variance and CI on the optimal cut-points are calculated but not
     returned in the final DiagTest3Grp object output.}     
   \item{partialDeriv}{A numeric data frame with one row and multiple
     columns, containing estimated SD parameters in each diagnosis group and the partial derivatives of Youden estimate
     w.r.t the relevant mean and SD parameters which are outputted for performance
     of statistical tests on markers under normal method or NA under
     nonparametric method.}
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
 \code{\link{VUS}}  \code{\link{Youden3Grp.PointEst}}  \code{\link{Youden3Grp.Variance.Bootstrap}} \code{\link{Youden3Grp.Variance.Normal}}   
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

 ##Normal method
 youden.res <- Youden3Grp(x=x, y=y, z=z, method = "Normal")
 print(youden.res)
 plot(youden.res)


 ##Kernel smoothing method with Sheather-Jones plug-in method
\dontrun{
 youden.res.KS <- Youden3Grp(x=x, y=y, z=z, method = "KS-SJ",randomStart.N=3)
    }
}

\keyword{htest}
\keyword{nonparametric}
\keyword{smooth}
\keyword{design}
\keyword{optimize}
