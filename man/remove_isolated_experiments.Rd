\name{remove_isolated_experiments}
\alias{remove_isolated_experiments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
remove_isolated_experiments(experiments, biological.group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experiments}{A list of wrapped experiments.}
  \item{biological.group}{A character string indicating the biological covariate that makes the biological groups. This must be the name of a column of the experiment.}
}
\details{
The experiments of the list that have no biological group in common with the other experiments are removed. This is an essential step to the integration of experiments using batch effect correction. Without it, the biological groups that are unshared between the experiments are confounding factors with the batch effect, and the latter cannot be corrected. The ultimate condition for this is the graph of biological intersections between the experiments to be a connected graph.

WARNING : this function only checks and removes the experiments that are isolated. However, there can remain several clusters of experiments that are part of a connected subgraph, while the graph of all the experiments isn't connected. The graph of intersections after the removal of experiments is displayed so that the user can notice it. In this condition, the user has to choose manually some experiments that make a connected subgraph.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
  \code{\link{merge_experiments}},
  \code{\link{correct_batch_effect}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
