\name{eigenangles}
\alias{eigenangles}
\title{Eigengenes Angles}
\usage{
eigenangles(data,batch,group,ref=NULL)
}
\description{
Calculates angles that quantifies the integration and the conservation of batches in an integrated experiment after batch effect correction.
}
\arguments{
  \item{data}{An expression matrix where genes correspond to rows and samples correspond to columns.}
  \item{batch}{A vector summarising the batch covariate, i.e. artifact factor.}
  \item{group}{A vector summarising biological groups that should explain differences between the distributions of the batches. To calculate the eigenangles between two batches, PCA is performed only on samples with biological properties common to both batches. }
  \item{reference}{An expression matrix with the same dimensions as \code{data} which corresponds to the reference dataset for the calculation of batch conservation angles. Typically, \code{data} corresponds to the corrected dataset and \code{ref} corresponds to the uncorrected dataset.}
}
\value{
A tibble with three columns :
  \item{$batch_}{This is a vector with the names of the batches.}
  \item{$integration_angles}{This is a list with the vectors of batch integration angles (each vector contains angles between the eigengenes of every rank of PCA).}
  \item{$conservation_angles}{This is a list with the vectors of batch conservation angles.}
}
\examples{
  
}
\author{Guillaume Heger, EMBL-EBI Trainee}
\references{\link{https://www.researchgate.net/publication/336533527_A_geometrical_approach_based_on_PCA_to_benchmark_the_algorithms_of_batch_effect_correction_applied_to_the_integration_of_RNA-Seq_data}}
