\name{eigenangles_wrapper}
\alias{eigenangles_wrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Eigengenes Angles for a wrapped integrated dataset}
\description{
Apply \code{eigenangles} to a wrapped integrated experiment.
}
\usage{
eigenangles_wrapper(experiment, biological.group, ref)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experiment}{A \code{SummarizedExperiment} object corresponding to a corrected integrated dataset.}
  \item{biological.group}{A character string giving the name of the column of \code{colData(experiment)} to use as biological covariate.}
  \item{reference}{A \code{SummarizedExperiment} corresponding to the uncorrected dataset.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A tibble of Eigengenes Angles.
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
\references{\url{https://www.researchgate.net/publication/336533527_A_geometrical_approach_based_on_PCA_to_benchmark_the_algorithms_of_batch_effect_correction_applied_to_the_integration_of_RNA-Seq_data}}
