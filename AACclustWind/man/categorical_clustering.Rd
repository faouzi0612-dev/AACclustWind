\name{categorical_clustering}
\alias{categorical_clustering}
\title{Categorical Clustering - Unified Interface}
\description{
  Run any of the supported categorical clustering methods and evaluate
  the partition using the Quality Index (QI).
}
\usage{
categorical_clustering(data, method = "aac_rcpp", k_range = NULL,
                       n_cores = NULL, seed = 42, verbose = FALSE, ...)
}
\arguments{
  \item{data}{A data.frame or matrix with categorical variables.}
  \item{method}{Character: one of "aac_sequential", "aac_parallel",
    "aac_rcpp", "hierarchical", "diana", "dbscan", "hdbscan",
    "spectral", "pop".}
  \item{k_range}{Range of k values to evaluate.}
  \item{n_cores}{Number of cores for parallel methods.}
  \item{seed}{Random seed.}
  \item{verbose}{Print progress.}
  \item{...}{Additional arguments.}
}
\value{
  A list with components: method, k_opt, QI, groups, time.
}
\examples{
\dontrun{
  data <- data.frame(a = sample(1:3, 20, TRUE),
                     b = sample(1:2, 20, TRUE))
  res <- categorical_clustering(data, method = "hierarchical")
}
}