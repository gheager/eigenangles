\name{merge_experiments}
\alias{merge_experiments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rough merger of experiments
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
merge_experiments(experiments, log, filter.unexpressed.genes, force=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experiments}{A list of wrapped experiments.}
  \item{log}{A logical indicationg whether the data shall be log-transformed.}
  \item{filter.unexpressed.genes}{A logical indicating if the genes that are unexpressed across all the samples of a batch shall be removed. This is for the proper functioning of the current version of ComBat in \code{correct_batch_effect}.}
  \item{force}{A logical. If the experiments of the list are of different classes, the function is stopped, unless \code{force=TRUE}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A wrapped experiment.
}
\references{

}
\author{Guillaume Heger, EMBL-EBI Trainee}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{load_experiments}},
  \code{\link{remove_isolated_experiments}},
  \code{\link{correct_batch_effect}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
