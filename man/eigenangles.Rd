\name{eigenangles}
\alias{eigenangles}
\title{Eigengenes Angles}
\usage{
eigenangles(data,batch,group,scale=FALSE,average=TRUE,verbose=TRUE)
}
\description{
Calculates angles between principal components of batches.
}
\arguments{
  \item{data}{A measure matrix where genes correspond to rows and samples correspond to columns.}
  \item{batch}{A vector summarising the batches covariate, i.e. artifact factor.}
  \item{group}{A vector summarising biological groups that should explain differences between the distributions of the batches. To calculate the eigenangles between two batches, PCA is performed only on samples with biological properties common to both batches. }
  \item{scale}{Logical. Argument for the \code{prcomp} function, determining whether the data shall be scaled before PCA. If \code{FALSE} (default), PCA analyses variance. If \code{TRUE}, PCA analyses correlation.}
  \item{average}{Logical. If \code{TRUE} (default and highly recommended in most cases), within each batch, samples from the same biological group are averaged in order to have the same number of samples from each batch so that the distributions can be assumed to be equals. Moreover it makes somehow that PCA analyses only the variance existing between the biological groups (which is the wanted information) and not variance coming from experimental processes.}
  \item{verbose}{Logical. Whether to display achieved comparisons.}
}
\value{
A list with two slots:
  \item{$batch_vs_all}{A list containing a vector of angles for each batch. Each vector contains the angles between the PCA subspaces for the considered batch and for the integrated dataset. There are as many angles as different biological groups in the batch if \code{average=TRUE}, and as many angles as samples in the batch if \code{average=FALSE}. Although there is no correspondance between the biological groups and the angles.}
  \item{$inter_batch}{A list of lists of vectors summarising all the comparisons between batches. Each vector contains the angles between the PCA subspaces for both batches in the comparison. There are as many angles as common biological groups between the batch if \code{average=TRUE}, and as many angles as samples in the smallest batch of the comparison if \code{average=FALSE}. Although there is no correspondance between the biological groups and the angles. If there is no common point between the batches, the vector is empty.}
}
\examples{

}
\author{Guillaume Heger, EMBL-EBI Trainee}
