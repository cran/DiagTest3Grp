\name{Pairwise.DiagTest3Grp.Test}
\alias{Pairwise.DiagTest3Grp.Test}

\title{
 Do pairwise statistical test on the summary measure (VUS or Youden)
 estimate between a pair among multiple markers 
}
\description{
 Do pairwise statistical test to compare the summary measure (VUS or Youden)
  estimates between a pair of markers among multiple markers assuming
  that the marker follows normal distribution in each diagnosis group, for both
  independent samples and paired samples.Resulting p-values can be
  adjusted for mulltiple testing. The function returns a nice
  upper-triagnular graphical display and print
  matrix on pairwise comparisons.
}

\usage{
Pairwise.DiagTest3Grp.Test(dat,paired=FALSE,type=c("VUS","Youden"),p=0,q=0,mu=0,
    conf.level=0.95, alternative=c("two.sided","less","greater"),
    p.adjust.method=c("none","holm","hochberg","hommel","bonferroni","BH","BY","fdr"),
    digits=3,...)
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
  \item{p.adjust.method}{A vector string specifying multiple testing
    adjustment method.Default=\dQuote{none}.}
  \item{digits}{A numeric integer. The number of digits for p-value
    print}
  \item{...}{Other arguments that can be input into heatmap.2() function
  in R package gplots.}
}

\details{
  see \emph{p.adjust} for details on multiple testing adjustment on p-values.
  and see \code{\link{DiagTest3Grp.Test}} for statistical test between
  two markers.
}
\value{
  Besides a heatmap illustrating the mulitple testing adjustment
  corrected P-values from pairwise comparisons in a upper-triangular
  display, the function also returns a list of two components:
  \item{print.matrix}{A data frame. If K markers are analyzed, the data
    frame will have K-1 row and K-1 columns. The first column contains the IDs
    of markers 1:(K-1) while the column names are ID for marker 2:K.The
    entries in column 3 and following are stat (t statistic), raw P
    (raw P value) or adjusted P  from testing two markers. The 2nd column indicates whether the
    matching row indicate stat (z-statistic), raw P
    (raw P value) or adjusted P. So this output basically presen the
    results in an upper triangular matrix format.}
  \item{pval.matrix}{More detailed output. A data frame with the
    columns: index of marker 1, index of marker 2, ID of marker1,summary
    measure estimated for marker 1, ID of
    marker2,summary
    measure estimated for marker 1, the lower limit of the resulting CI
    on the difference between the summary measure estimates (marker1-marker2), the upper
    limit of the CI, z-statistic, raw P value and adjusted P value.}
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

  Luo, J and Xiong, C. (2012) Youden Index and Associated Optimal Cut-point
  for Three Ordinal Groups. \emph{Communications In Statistics-Simulation
  and Computation}(in press).

  Jingqin Luo, Chengjie Xiong (2012). DiagTest3Grp: An R Package for
  Analyzing Diagnostic Tests with Three Ordinal Groups.\emph{Journal of
    Statistical Software}, \bold{51}(3), 1-24.URL
  http://www.jstatsoft.org/v51/i03/.
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
 new.AL <- data.frame(group=AL$group,-AL[,-1])
 new.AL <- new.AL[,1:4]##use the first 3 markers as an example

 ##omnibus test on the 3 AL markers
 Pairwise.DiagTest3Grp.Test(dat=new.AL, paired=TRUE, type="VUS",p=0,q=0,
               conf.level=0.95,p.adjust.method="fdr",digits=3)

}
\keyword{htest}
