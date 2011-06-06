\name{AL}
\alias{AL}
\docType{data}
\title{Alzheimer's disease neuropsychometric marker dataset
}
\description{
  The dataset is a subset of the longitudinal cohort of Washington University (WU)
  Alzheimer's Disease Research Center (ADRC). In the AL dataset, measurements of 14 neuropsychological markers
  were collected on 118 independent individuals of age 75. The individuals were classified into 3 groups based on
  publsihed clinical demential rating (CDR).   
}
\usage{data(AL)}
\format{
  A data frame with 118 observations on the following 15 variables.
  \describe{
    \item{\code{group}}{a factor with levels \eqn{D^-}, \eqn{D^0}
      \eqn{D^+}, the \eqn{D^-} group has CDR (clinical demential
      rating)=0, i.e., no dementia group; the \eqn{D^0} group has
      CDR=0.5, i.e., mild dementia and the \eqn{D^+} group has CDR=1,
      Alzheimer's disease. 
    }
    \item{\code{FACTOR1}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{global factor}.} 
    \item{\code{ktemp}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{temporal factor}.}
    \item{\code{kpar}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{parietal factor}.}
    \item{\code{kfront}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{frontal factor}.} 
    \item{\code{zpsy004}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{longical memory}.}
    \item{\code{zpsy005}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{digital span forward}.}
    \item{\code{zpsy006}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{digital span backward}.} 
    \item{\code{zinfo}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{information}.} 
    \item{\code{zbentc}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{visual retention (10s)}.}
    \item{\code{zbentd}}{a numeric vector, measurements on the
      neuropsychometric test for \dQuote{visual retention (copy)}.} 
    \item{\code{zboston}}{a numeric vector,a numeric vector,
      measurements on the neuropsychometric test for \dQuote{boston
      naming}.} 
    \item{\code{zmentcon}}{a numeric vector,measurements on the
      neuropsychometric test for \dQuote{mental control}.}
    \item{\code{zworflu}}{a numeric vector,measurements on the
      neuropsychometric test for \dQuote{word fluency}.}
    \item{\code{zassc}}{a numeric vector,measurements on the
      neuropsychometric test for \dQuote{associate learning}.}
  }
}
\details{
 See references for more details on the explanation of the 14 markers
 and scoring details.  
}

\source{
  Xiong, C. and van Belle, G. and Miller, J.P. and Morris, J.C. (2006)
  Measuring and Estimating Diagnostic Accuracy When There 
  Are Three Ordinal Diagnostic Groups. \emph{Statistics In Medicine} \bold{25}
  7 1251--1273.
  
  Sorandt M, Hill RD. (1989) Very mild snile dementia of the Alzheimer type. II
  Psychometric test performance. \emph{Archives of Neurology} \bold{46} 383--386. 
}
\examples{
    data(AL)

    ## maybe str(AL) ; plot(AL);  ...
}
\keyword{datasets}
