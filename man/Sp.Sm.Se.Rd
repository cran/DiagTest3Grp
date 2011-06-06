\name{Sp.Sm.Se}
\alias{Sp.Sm.Se}

\title{
Obtain empirical estimation of the three correct classification probabilities
}
\description{ once identified the optimal cut-points from VUS or
  Youden3Grp analyses, use the function to obtain the three empirical correct
  classification probabilities associated with each diagnosis group. 
}

\usage{
  Sp.Sm.Se(x, y, z, t.minus, t.plus)
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
 \item{t.minus}{
 A numeric value, the lower optimal cut-point identified from VUS or Youden3Grp
   analyses.
 }
 \item{t.plus}{
 A numeric value,  the upper optimal cut-point identified from VUS or
   Youden3Grp analyses.
  }
 }
 \details{
   Specificity:  \eqn{Sp=Pr(x\le t_-)} for \eqn{D^-} group;
   Sensitivity: \eqn{Se=Pr(z\ge t_+)} for \eqn{D^+} group and the probability of the diagnostic
   test for the \eqn{D^0} group fall between the two cut points:
   \eqn{Sm=Pr(t_- \le y \le t_+)}.
   These three probabilities will be estimated empirically.
   
 }
\value{
  Return a numeric vector with three components Sp, Sm and Se, the three
  correct classification probabilities.
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
 \code{\link{VUS}}  \code{\link{Normal.VUS}} \code{\link{NonParametric.VUS.var}}
}
\examples{

 data(AL)
 group <- AL$group

 ##take the negated kfront marker measurements
 kfront <- -AL$kfront

 x <- kfront[group=="D-"]
 y <- kfront[group=="D0"]
 z <- kfront[group=="D+"]
 class.prob <- Sp.Sm.Se(x=x,y=y,z=z,t.minus=-1.68,t.plus=0.91)
}


