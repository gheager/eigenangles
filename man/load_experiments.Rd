\name{load_experiments}
\alias{load_experiments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
load_experiments(directory, item.SimpleList='rnaseq')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{directory}{A character string that is the path to the directory where all the experiments to load, and not any other file, are stored in Rdata format.}
  \item{item.SimpleList}{If the experiments are wrapped in a \code{SimpleList} object, this is the accessor (character string or integer) of the actual experiment object within this \code{SimpleList}. The experiment is unwrapped from this \code{SimpleList}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of the loaded experiments.
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
  \code{\link{remove_isolated_experiments}},
  \code{\link{merge_experiments}},
  \code{\link{correct_batch_effect}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
