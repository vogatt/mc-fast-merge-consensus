\name{mcmergecons}
\alias{mcmergecons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multicore enabled version of merged consensus clustering.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function performs merged consensus clustering.
It allows for usage of multiple clustering algoritms.
Set number of cores with the numofcores parameter.
}
\usage{
mcmergecons(x, diss = FALSE, algorithms = list("agnes"), alparams = list(), alweights = list(), clmin = 2, clmax = 10, prop = 0.8, reps = 50, merge = 0, numofcores = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Data frame containing data or data frame with distance matrix.
%%     ~~Describe \code{x} here~~
}
  \item{diss}{
  TRUE if x is a distance matrix.
  FALSE otherwise.
%%     ~~Describe \code{diss} here~~
}
  \item{algorithms}{
  List containing inner algorithms names.
%%     ~~Describe \code{algorithms} here~~
}
  \item{alparams}{
  List passing parameters to inner algorithms.
%%     ~~Describe \code{alparams} here~~
}
  \item{alweights}{
  List passing weights for merging.
%%     ~~Describe \code{alweights} here~~
}
  \item{clmin}{
  Minimum number of clusters expected in the dataset.
%%     ~~Describe \code{clmin} here~~
}
  \item{clmax}{
  Maximum number of clusters expected in the dataset.
%%     ~~Describe \code{clmax} here~~
}
  \item{prop}{
  Percentage of data used in subsampling. Must be in [0.0 , 1.0] range.
%%     ~~Describe \code{prop} here~~
}
  \item{reps}{
  Number of repeted runs of algorithms to perform.
%%     ~~Describe \code{reps} here~~
}
  \item{merge}{
  If more than 0.0 merging is performed.
%%     ~~Describe \code{merge} here~~
}
  \item{numofcores}{
  Number of available cores.
%%     ~~Describe \code{numofcores} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list with results for each expected number of clusters and algorithm.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
T Ian Simpson, J Douglas Armstrong, and Andrew Jarman. Merged consensus clustering
to assess and improve class discovery with microarray data. BMC bioinformatics, 11(1):590,
2010.
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
a <- mcmergecons(x=data.frame(data),diss=FALSE,algorithms=list('kmeans','agnes'),numofcores=2,merge=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
