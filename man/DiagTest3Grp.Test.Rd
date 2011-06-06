\name{DiagTest3Grp.Test}
\alias{DiagTest3Grp.Test}

\title{
 Statistical test on the summary measure (VUS or Youden) estimate on two
 markers or multiple markers
}
\description{
  Statistical test to compare the summary measure (VUS or Youden)
  estimates between two markers or multiple markers assuming that the marker
  follows normal distribution in each diagnosis group, for both
  independent samples and paired samples.
}

\usage{
DiagTest3Grp.Test(dat,paired=FALSE,type=c("VUS","Youden"),p=0,q=0,mu=0,conf.level=0.95,alternative=c("two.sided","less","greater"))
}

\arguments{
  \item{dat}{For unpaired data (data measured on different independent
    sets of samples), dat is a list of K (=number of markers)
    components, each corresponding to a marker. Each component of the 
    list should be a data frame of number of observation by two columns:(1) group
    (note:must be at the 1st column) with the three levels of \eqn{D^-},
    \eqn{D^0} and \eqn{D^+} (corresponding to increasing disease
    severity respectively);  (2)numeric marker measurements at the 2nd 
    column.
    for paired data, dat should be a data frame with
    observations at row and group/markers at column (group must be at the
    1st column) with the three levels \eqn{D^-},
    \eqn{D^0} and \eqn{D^+}(corresponding to increasing disease
    severity respectively), followed by markers at the 2nd to the
    last column.}
  \item{paired}{A logic value TRUE/FALSE, indicating whether all markers are measured on the same set (paired=TRUE) or not (paired=FALSE).}
  \item{type}{A character value. type=\dQuote{VUS} or \dQuote{Youden},
    indicates which summary measure (VUS or Youden) will be used to
    indicate the discriminative ability of the three-group diagnostic
    tests. Allows unique argument matching,i.e., can use unique
    abbreviations such as \dQuote{V} or \dQuote{Y}.}
  \item{p}{A numeric value. To test on partial VUS satisfying specificity >=p and
    sensitivity >=q, set p, q to desired value. default=0 for full VUS.}
  \item{q}{A numeric value.To test on partial VUS satisfying specificity
    >=p and sensitivity >=q, set p, q to desired value. default=0 for
    full VUS.}
  \item{mu}{A numeric value. The hypothesized true difference between
    two markers. Default mu=0 but can be changed to a nonzero value.}
  \item{conf.level}{a numerica value. The confidence
    interval. default=95\%. The confidence interval on point estimate of
    abs diff between two markers's summary measures will be
    returned. For >=3 markers, this will be NA.}
  \item{alternative}{A character value. Now  is fixed as "two-sided".}   
}
\details{
  Let \eqn{\theta_j} indicate the summary measure (either VUS or Youden)
  on a marker j. When testing Q markers (\eqn{Q\ge3}) markers, the omnibus test similar to one-way ANOVA is
  performed with the hypothesis {H_0}:
  \eqn{\theta_1=\theta_2=\cdots=\theta_Q} against the
  alternative \eqn{H_1}: at least a pair of markers have unequivalent
  summary measures. The resulting statistic is a Chi-square statistic
  with Q-1 degree of freedom. When testing 2 markers, the normal test is
  performed to test the hypothesis {H_0}:
  \eqn{\theta_1-\theta_2=\theta_0} against
\eqn{H_1}:\eqn{\theta_1-\theta_2\neq 0} for
\eqn{alternative=\dQuote{two.sided}},
\eqn{H_1}:\eqn{\theta_1-\theta_2<0} for \eqn{alternative=\dQuote{less}}
and \eqn{H_1}:\eqn{\theta_1-\theta_2> 0} for \eqn{alternative=\dQuote{greater}}  The z-statistic is
  formed (difference in summary measures divided by the standard error
  on the difference) for normal test.
}
\value{
  An object of the classic S3 class \emph{htest} (slightly modified) with the list of
  components.
  \item{statistic}{A numeric value. In situations of testing >=3
    markers, this is a Chi-Square statistic; in testing 2 markers, this
    is a Z-statistic. }
  \item{p.value}{A numeric value. The p-value testing the hypothesis
    (see details) from Chi-square test (for >=3 markers) or normal test
    (2 markers).}
  \item{parameter}{A numeric value. For the Chi-square test on >=3
    markers, it's the degree of freedom. In the normal test on 2
    markers, it's mean 0 and sd=1.}
  \item{estimate}{A numeric vector. The summary measures estimated on
    each marker.}
  \item{conf.int}{A numeric vector of length 2. The confidence interval
    on the difference in the summary measures of two marker(marker1 at
    the 1st column of dat input-marker2 at the 2nd column at dat
    input). For testing >=3 markers, NA.}
  \item{method}{A character value. \dQuote{Chi-square test} for >=3 markers or \dQuote{Normal-test}
    for 2 markers.}
  \item{type}{A character value. Either VUS or Youden, which type of
    summary measure is in use.}
  \item{null.value}{A numeric value. The inputted value for mu which is now fixed as 0 for >=3 marker
    comparison but can be changed to any numeric value for two marker comparison.}
  \item{alternative}{A character value which is always \dQuote{two.sided} for
    >=3 marker comparison but can be \dQuote{two.sided}, \dQuote{less}
    or \dQuote{greater}.}
  \item{data.name}{Usually prints Test of the name of the inputted data
    frame on the type of a summary measure.}
  \item{Sigma}{A numeric matrix. The covariance matrix on the summary
    measures estimated on markers.}
   }
  
\references{
   Xiong, C. and van Belle, G. and Miller, J.P. and Morris, J.C. (2006)
  Measuring and Estimating Diagnostic Accuracy When There 
  Are Three Ordinal Diagnostic Groups. \emph{Statistics In Medicine} \bold{25}
  7 1251--1273.

  Ferri, C. and Hernandez-Orallo, J. and Salido, M.A. (2003) Volume
  under the ROC Surface for Multi-class Problems \emph{LECTURE NOTES IN
    COMPUTER SCIENCE} 108--120. 
  Xiong, C and Van Belle, G. and Miller, J.P. and Yan, Y. and Gao, F.,
  and Yu, K, and Morris, J. (2007) A Parametric Comparison of Diagnostic
  Accuracy with Three Ordinal Diagnostic Groups. \emph{Biometrical
    Journal} 5, 682-693

  Luo, J and Xiong, C. (2011) Youden Index and Associated Optimal
  Cut-points for Three Ordinal Groups.\emph{Communications in
    Statistics--simulation and computation, submitted}
  
  Luo, J and Xiong, C. (2011) DiagTest3Grp: An R package for Analyzing
  Diagnostic Tests with Three Ordinal Groups.\emph{Journal of
    Statistical Software, submitted}  
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
 \code{\link{VUS}} \code{\link{Youden3Grp}}  
}
\examples{

 data(AL)
 AL <- data.frame(group=AL$group,-AL[,-1])

 #group <- AL$group
 #table(group)

 ##omnibus test on all 14 marker
 DiagTest3Grp.Test(dat=AL,paired=TRUE,type="VUS",p=0,q=0,conf.level=0.95)
 DiagTest3Grp.Test(dat=AL,paired=TRUE,type="Youden",p=0,q=0,conf.level=0.95)

 ##normal test on 2 of the markers
 new.AL <- subset(AL,select=c(group,ktemp,zbentd))
 DiagTest3Grp.Test(dat=new.AL,paired=TRUE,type="VUS",p=0,q=0,conf.level=0.95)
 DiagTest3Grp.Test(dat=new.AL,paired=TRUE,type="Youden",p=0,q=0,conf.level=0.95)

}
\keyword{htest}
