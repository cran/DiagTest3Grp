\name{Youden3Grp.PointEst}
\alias{Youden3Grp.PointEst}
 
\title{
  Provide point estimates for the extended Youden index, the associated
  optimal cut-points, and sample size for three ordinal groups
}
\description{
 Estimate for the extended Youden index, the associated optimal
 cut-points using normal and non-parametric methods and the
 associated. The sample size  to estimate the extended Youden index of a
 diagnostic marker within a given margin of error under normal
 assumptions is also provided.   
}
\usage{ Youden3Grp.PointEst(x, y, z, method = "Normal", randomStart.N
=1, optim.method = NULL, t.minus.start = NULL, t.plus.start = NULL, 
lam.minus = 1/3, lam0 = 1/3, lam.plus = 1/3, typeIerror = 0.05, margin=0.05,FisherZ=FALSE, ...)  
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
    cut-point (t-) which separate the \eqn{D^-} and \eqn{D^0} group. Default
    t.minus.start=NULL will randomly generate  starting value. }
\item{t.plus.start}{A numeric value. The starting points of the upper optimal
    cut-point (t+) which separate the \eqn{D^0} and \eqn{D^+}
    group. Default t.plus.start=NULL will randomly generate starting value. }
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
  \item{FisherZ}{Reference to the argument in \code{\link{Normal.VUS}}.}
 \item{...}{ Other options for the R function \code{\link{optim}}.}
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
  a list of three components. \dQuote{dat}: a list of compoenents x, y, z,
  recording raw data inputted (with NAs removed); \dQuote{dat.summary}:
  a data frame of three rows (corresponding to \eqn{D^-}, \eqn{D^0} and
  \eqn{D^+}) separately and three columns (number of non-missing observations, mean
  and SD of raw data); \dQuote{est}:  a data frame of one row and the following
  columns numeric value components:
  
  \item{mu.minus}{The sample mean estimate of the a diagnostic's
    measurements in the \eqn{D^-} group, which is sample mean of data if
    method=\dQuote{Normal}/\dQuote{EMP}/\dQuote{KS}/\dQuote{KS-SJ} while
    if method=\dQuote{TN}, mu.minus is the sample mean of 
    the data after Box-cox transformation. }
  \item{mu0}{The sample mean estimate of the a diagnostic's
    measurements in the D0 group, which is sample mean of data if
    method=Normal while if method=TN, mu0 is the sample mean of
    the data after Box-cox transformation. }
  \item{mu.plus}{The sample mean estimate of the a diagnostic's
    measurements in the D+ group, which is sample mean of data if
    method=Normal while if method=TN, mu.plus is the sample mean of
    the data after Box-cox transformation. }
  \item{s.minus}{The sample standard deviation (SD) estimate of the a diagnostic's
    measurements in the D- group, which is sample SD of data if
    method=Normal while if method=TN, the s.minus is the sample SD of
    the data after Box-cox transformation. }
  \item{s0}{The sample standard deviation (SD) estimate of the a diagnostic's
    measurements in the D0 group, which is sample SD of data if
    method=Normal while if method=TN, the s0 is the sample SD of
    the data after Box-cox transformation. }
  \item{s.plus}{The sample standard deviation (SD) estimate of the a diagnostic's
    measurements in the D+ group, which is sample SD of data if
    method=Normal while if method=TN, the s.plus is the sample SD of
    the data after Box-cox transformation. }
  \item{t.minus}{The estimate of the lower optimal cut-point.}
  \item{t.plus}{The estimate of the upper optimal cut-point.}
  \item{Se}{Sensitivity, the probability of a randomly selected subject
    from D+ group has a test result above t.plus.  }
  \item{Sp}{Specificity the probability of a randomly selected subject
    from D- group has a test result below t.minus.  }
  \item{Sm}{The probability of a randomly selected subject
    from D0 group has a test result between t.minus and t.plus.  }
 \item{youden}{The extended Youden index for three ordinal groups
   attained at the optimal cut-point t.minus and t.plus.  }
  \item{youden.z}{The Fisher's Z transformation of the extended Youden
    index.  }
  \item{lambda.est}{The estimated lambda, the parameter used in
    Box-Cox transformation y=(x^lambda-1)/lambda. The estimate will be a
    numeric value when the method=TN is
    used to estimate the extended Youden index and NA when other methods
    were used.  }
  \item{t.minus.TN}{The optimal lower cut-point t.minus estimated for the data after Box-Cox
    transformation. The estimate will be a numeric value when the method=TN is
    used to estimate the extended Youden index and NA when other methods
    are used. }   
  \item{t.plus.TN}{The optimal upper cut-point t.plus estimated for the
    data after Box-Cox transformation. The estimate will be a numeric
    value when the method=TN is used to estimate the extended Youden
    index and NA when other methods are used.} 
  \item{sampleSize}{The sample size on investigating a diagnostic test
     based on pilot measurments x, y, z. The sample size calculated such
     that based on given population proportion of the three diagnosis
     groups, the margin of error on the (1-typeIerror)*100\% confidence
     interval of the extended Youden index estimate under normality is
     as specified, by default 0.05. }   
}

\references{
  Luo, J and Xiong, C. (2012) Youden Index and Associated Optimal Cut-point
  for Three Ordinal Groups. \emph{Communications In Statistics-Simulation
  and Computation} (in press).
}
\author{
 Jingqin Luo
}
\note{
 Bug reports, malfunctioning, or suggestions for further improvements or
 contributions can be sent to Jingqin Luo <rosy@wubios.wustl.edu>.
}


\seealso{
 \code{\link{Youden3Grp.Variance.Bootstrap}} 
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
 youden.res.normal <- Youden3Grp.PointEst(x=x, y=y, z=z, method = "Normal")
 youden.res.normal$est$youden

 ##Kernel smoothing method with Sheather-Jones plug-in method
 youden.res.KS <- Youden3Grp.PointEst(x=x, y=y, z=z, method = "KS-SJ",randomStart.N=3)
 youden.res.KS$est$youden    
}

\keyword{htest}
\keyword{nonparametric}
\keyword{smooth}
\keyword{design}
\keyword{optimize}
