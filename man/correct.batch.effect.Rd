\name{correct.batch.effect}
\alias{correct.batch.effect}
\title{
Batch effect correction using several methods
}
\description{
Corrects batch effect within an integrated dataset.}
\usage{
correct.batch.effect(data, batch, method=c('none','combat','ruv','mnn'), 
                     model, model.data, log=TRUE, k=1)
}
\arguments{
  \item{data}{A count matrix with genes in rows and samples in columns.}
  \item{batch}{A vector indicating batches of samples.}
  \item{method}{A string indicating which method to use for batch effect correction. Choices are \code{"none"} for no correction, \code{"combat"} for the \code{ComBat} function from package \code{sva}, \code{"ruv"} for \code{RUVs} function from package \code{RUVSeq} and \code{"mnn"} for \code{mnnCorrect} function from packge \code{batchelor}.}
  \item{model}{A formula indicating a biological model to take into account during batch effect correction.}
  \item{model.data}{A model frame containing data for the \code{model} argument.}
  \item{log}{A logical indicating whether to log transform the count matrix before correction. Default is \code{TRUE}. NB: if \code{method="ruv"}, data are log transformed anyway.}
  \item{k}{An integer parameter used only for \code{method="ruv"} or \code{method="mnn"}. See \code{\link[RUVSeq]{RUVs}} and \code{\link[batchelor]{mnnCorrect}} for more information.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A corrected expression matrix. NB: depending on the value of \code{log} argument, the elements of the matrix should be considered either as log-counts or as counts.
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
  \code{\link{integrate.experiments}}
  \code{\link[sva]{ComBat}}
  \code{\link[RUVSeq]{RUVs}}
  \code{\link[batchelor]{mnnCorrect}}
}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
