\name{integrate.experiments}
\alias{integrate.experiments}
\title{
Experiments integration with batch effect correction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
integrate.experiments(...,list=NULL,model=NULL,
                     method=c('none','combat','ruv','mnn'),log=TRUE,k=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{As many \code{SummarizedExperiment} objects as wanted. Every argument passed into \code{...} will be considered as a batch.}
  \item{list}{An alternative to \code{...} if the \code{SummarizedExperiment} objects are already gathered in a list. Every element of the list will be considered as a batch.}
  \item{model}{A formula indicating a biological model to take into account during batch effect correction. The terms of the formula must be columns names in the \code{colData()} field of every experiment given as argument.}
  \item{method}{A string indicating which method to use for batch effect correction. Choices are \code{"none"} for no correction, \code{"combat"} for the \code{ComBat} function from package \code{sva}, \code{"ruv"} for \code{RUVs} function from package \code{RUVSeq} and \code{"mnn"} for \code{mnnCorrect} function from packge \code{batchelor}.}
  \item{log}{A logical indicating whether to log transform the count matrix before correction. Default is \code{TRUE}. NB: if \code{method="ruv"}, data are log transformed anyway.}
  \item{k}{An integer parameter used only for \code{method="ruv"} or \code{method="mnn"}. Type \code{?sva::ComBat} and \code{?batchelor::mnnCorrect} in R console for more information.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A \code{SummarizedExperiment} object whose \code{assays()} field contains a \code{$corrected} element which is the corrected expression matrix. NB: depending on the value of \code{log} argument, the elements of the matrix should be considered either as log-counts or as counts.

Its \code{colData()} field is the merger of all the \code{colData()} of the experiments given as arguments.

Its \code{metadata()} field contains a \code{$batch} element giving the batches for every sample.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
