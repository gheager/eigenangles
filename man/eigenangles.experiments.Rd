\name{eigenangles.experiments}
\alias{eigenangles.experiments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Eigengenes angles between the batches of an integrated experiment}
\description{
Evaluates batch effect correction by using eigengenes angles.}
\usage{
eigenangles.experiments(experiment,group,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experiment}{A \code{SummarizedExperiment} object resulting from a call to \code{\link{integrate.experiments}} function.}
  \item{group}{A string indicating a biological group to take into account to compare batches wisely. Must be a column name of the \code{colData()} field of the \code{experiment} argument.}
  \item{...}{Any other "technical" argument to be provided to \code{\link{eigenangles}} function.}
}
\details{
  \code{eigenangles.experiments(experiment,group,...)} is equivalent to \code{eigenangles(data=assays(experiment)$corrected, batch=metadata(experiment)$batch, group=experiment[[group]], ...)}
}
\value{
  The output of a call to \code{\link{eigenangles}} function.
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
  \code{\link{integrate.experiments}}
}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
