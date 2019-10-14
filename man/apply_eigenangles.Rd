\name{apply_eigenangles}
\alias{apply_eigenangles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Eigengenes Angles to benchmark methods of batch effect correction
}
\description{
Apply \code{eigenangles} to a set of corrected datasets from different correction methods.
}
\usage{
apply_eigenangles(..., uncorrected, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{Multiple arguments corresponding to as many corrected datasets (\code{SummarizedExperiment} objects) as wanted to be compared. These arguments should be named according to the correction method used, in order to identify them easily in the result.}
  \item{uncorrected}{A dataset (\code{SummarizedExperiment}) that is the uncorrected dataset. This has to be declared independantly to be identified as the reference for angles of batch conservation.}
  \item{group}{A character string giving the name of the column of the \code{colData()} field of each experiment to use as biological covariate.}
}
\details{
Arguments in the ellipsis \code{...} can be \code{SummarizedExperiment} datasets or lists of \code{SummarizedExperiment} datasets. A list of \code{SummarizedExperiment} datasets can be used to gather datasets that have been corrected using the same method but with different values of a parameter k (if applicable to this correction method). This k parameter will be displayed in the \code{$k} column of the result, if applicable.
}
\value{
A tibble of four or five columns :
  \item{$algorithm}{The algorithm used to correct batch effect : they are taken from the names of the arguments in \code{...}. One of them refers to \code{uncorrected}.}
  \item{$k}{The value of the parameter k of the algorithm if applicable.}
  \item{$batch_}{This is a vector with the names of the batches.}
  \item{$integration_angles}{This is a list with the vectors of batch integration angles (each vector contains angles between the eigengenes of every rank of PCA).}
  \item{$conservation_angles}{This is a list with the vectors of batch conservation angles.}
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
  \code{\link{eigenangles}},
  \code{\link{do_eigenangles}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
