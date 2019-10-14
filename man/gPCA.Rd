\name{gPCA}
\alias{gPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Guided PCA}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gPCA(data,batch,scaleY=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An expression matrix with genes as rows and samples as columns.}
  \item{batch}{A vector indicating the groups to which the samples belong.}
  \item{scaleY}{A logical indicating if the batches should be averaged or summed. Default is \code{TRUE} which means the batches are averaged.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class `gpca`, which is a list containing the delta statistics, the variance ranks of the guided principal components, as well as the `prcomp` objects for guided and unguided PCA.
}
\references{
Sarah E. Reese, Kellie J. Archer, Terry M. Therneau, Elizabeth J. Atkinson, Celine M. Vachon, Mariza de Andrade, Jean-Pierre A. Kocher, Jeanette E. Eckel-Passow, A new statistic for identifying batch effects in high-throughput genomic data that uses guided principal component analysis, Bioinformatics, Volume 29, Issue 22, 15 November 2013, Pages 2877–2883, https://doi.org/10.1093/bioinformatics/btt480
}
\author{Guillaume Heger, EMBL-EBI Trainee}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
