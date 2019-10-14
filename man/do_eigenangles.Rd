\name{do_eigenangles}
\alias{do_eigenangles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Eigengenes Angles for a wrapped integrated dataset}
\description{
Apply \code{eigenangles} to a wrapped integrated experiment.
}
\usage{
do_eigenangles(experiment, group, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experiment}{A \code{SummarizedExperiment} object corresponding to a corrected integrated dataset.}
  \item{group}{A character string giving the name of the column of \code{colData(experiment)} to use as biological covariate.}
  \item{reference}{A \code{SummarizedExperiment} corresponding to the uncorrected dataset.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A tibble of Eigengenes Angles.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Guillaume Heger, EMBL-EBI Trainee}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{eigenangles}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
