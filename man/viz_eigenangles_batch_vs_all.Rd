\name{viz_eigenangles_batch_vs_all}
\alias{viz_eigenangles_batch_vs_all}
\title{
Visualising Eigengenes Angles between batches and integrated dataset
}
\description{
Plots eigengenes angles between each batch and integrated dataset.}
\usage{
viz_eigenangles_batch_vs_all(...)
}
\arguments{
  \item{...}{As many eigenangles lists as wanted for graphical comparison. All these lists must be results of \code{eigenangles} function. These lists shall have common batches. Giving a name to each argument may be useful to display names of the integrated datasets in the legend. Otherwise numbers are displayed, corresponding to the order of arguments in the function call.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A ggplot object. Displays as many plots as different batches in all the datasets.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Guillaume Heger, EMBL-EBI Trainee
}
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
