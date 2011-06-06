\name{VUS.CutPoint}
\alias{VUS.CutPoint}

\title{
 Find the optimal cut-points in volumn under ROC surface (VUS) analyses
}
\description{
   Find the optimal cut-points in volumn under ROC surface (VUS) analyses
 which correspond to the coordinate triplet with the shortest distance to the perfect classification
 coordinates (1,1,1).
}
\usage{
 VUS.CutPoint(x, y, z, cut.seq=NULL,optimize=F)
}

\arguments{
  \item{x}{
A numeric vector. A diagnostic test's measurements in  the D- (usually healthy subjects).
}
\item{y}{
A numeric vector.  A diagnostic test's measurements in  the D0 (usually mildly diseased subjects).
 }
 \item{z}{
A numeric vector.  A diagnostic test's measurements in  the D+ (usually severely diseased subjects).
 }
  \item{cut.seq}{ A numeric vector. A sequence of user-specified values from which
    the optimal cut-points will be chosen. If set to NULL, the unique
    collection of values of x,y,z will be used as the sequence.
  }
  \item{optimize}{A logical value of TRUE or FALSE. If FALSE, take the
    empirical optimal cut point identified by empirical search within the given
    cut.seq as final reported optimal cut point.If TRUE, using the
    empirical optimal cut point as starting point in optimization algorithm for final
    optimal cut point.}
  
}
\details{
   The optimal cut-point t.minus, t.plus is determined based on a pair
   of values \eqn{(t_-,t_+)}, \eqn{t_-<=t_+}, from cut.seq to minimize the goal function:
   \eqn{(Pr(x\le t_-)-1)^2+(Pr(t_-<y<t_+)-1)^2+(Pr(z\ge t_+)-1)^2}. The
    probabilities are estimated by empirical CDF. 
}
\value{
  Return a data frame with two columns \dQuote{t.minus} and
  \dQuote{t.plus} as optimal pair of cut-points for ROC surface analysis.
}

\author{
  Jingqin Luo
}
\note{
 Bug reports, malfunctioning, or suggestions for further improvements or
 contributions can be sent to Jingqin Luo <rosy@wubios.wustl.edu>.   
}

\seealso{
  \code{\link{VUS}}  \code{\link{Normal.VUS}} \code{\link{Youden3Grp}}
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

 vus.cut <- VUS.CutPoint(x=x,y=y,z=z,cut.seq=NULL)
}

\keyword{optimize}


