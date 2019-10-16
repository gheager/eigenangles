\name{correct_batch_effect}
\alias{correct_batch_effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Batch effect correction between roughly merged batches
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correct_batch_effect(experiment, model, method=c('ComBat','RUV','MNN'), k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experiment}{A wrapper of an experiment (\code{SummarizedExperiment} or \code{ExpressionSet} object).}
  \item{model}{A formula giving the biological model to take into account during batch effect correction. The possible variables are the columns names of the experiment.}
  \item{method}{A character string among 'ComBat', 'RUV' and 'MNN' indicating which correction method to use to correct batch effect.}
  \item{k}{An integer or a vector of integers that is the hyper-parameter for RUV and MNN. For RUV, this is the number of componants of unwanted variation to remove, while for MNN, it is the number of considered nearest neighbours per sample. If \code{k} is a vector of integers, the chosen method is run for each specified value of \code{k} and the result is the list of the individual results for each value of \code{k}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A wrapped experiment of the same format than \code{experiment} (or a list of such objects if \code{k} is a vector.)
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
  \code{\link{load_experiments}},
  \code{\link{remove_isolated_experiments}},
  \code{\link{merge_experiments}}

  \code{\link[sva]{ComBat}},
  \code{\link[RUVSeq]{RUVs}},
  \code{\link[batchelor]{mnnCorrect}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
