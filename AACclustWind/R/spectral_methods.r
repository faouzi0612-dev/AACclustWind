#' Spectral Clustering
#'
#' @param data A data.frame with categorical variables
#' @param k_range Range of k to evaluate for eigengap
#' @param seed Random seed
#' @return List with results
#' @export
spectral_clustering <- function(data, k_range = 2:10, seed = 42) {
  data_mat <- as.matrix(data)
  m <- nrow(data_mat)

  # Similarity matrix W
  W <- matrix(0, m, m)
  for (i in 1:m) {
    for (j in i:m) {
      val <- sum(data_mat[i, ] == data_mat[j, ])
      W[i, j] <- val
      W[j, i] <- val
    }
  }
  diag(W) <- 0

  # Eigengap
  D <- diag(rowSums(W))
  L <- D - W
  eigenvalues <- sort(base::eigen(L, symmetric = TRUE, only.values = TRUE)$values)

  gaps <- c()
  for (i in seq_along(k_range)) {
    if (k_range[i] < length(eigenvalues)) {
      gaps <- c(gaps, eigenvalues[k_range[i] + 1] - eigenvalues[k_range[i]])
    }
  }

  valid_gaps <- gaps[1:min(length(k_range) - 1, length(gaps))]
  k_opt <- k_range[which.max(valid_gaps)]

  # Normalized cut
  D_inv_sqrt <- diag(1 / sqrt(pmax(rowSums(W), 1e-10)))
  L_sym <- D_inv_sqrt %*% W %*% D_inv_sqrt
  eig <- base::eigen(L_sym, symmetric = TRUE)
  U <- eig$vectors[, 1:k_opt]
  row_norms <- sqrt(rowSums(U^2))
  U_norm <- U / pmax(row_norms, 1e-10)

  set.seed(seed)
  km <- stats::kmeans(U_norm, centers = k_opt, nstart = 20)
  grupos <- split(1:m, km$cluster)
  qi <- calculate_qi(data, grupos)

  list(
    method = "Spectral",
    k_opt = k_opt,
    QI = qi,
    groups = grupos,
    kmeans_obj = km
  )
}