#' Calculate Quality Index (QI) for a categorical partition
#'
#' @param data A data.frame or matrix with categorical variables
#' @param groups A list of integer vectors representing cluster memberships
#' @return Numeric QI value
#' @export
calculate_qi <- function(data, groups) {
  if (length(groups) == 0) return(NA)
  data_mat <- as.matrix(data)
  m <- nrow(data_mat)
  n <- ncol(data_mat)

  # Coincidence matrix C
  C <- matrix(0, m, m)
  for (i in 1:m) {
    for (j in i:m) {
      val <- sum(data_mat[i, ] == data_mat[j, ])
      C[i, j] <- val
      C[j, i] <- val
    }
  }

  # Membership matrix S
  S <- matrix(0, m, m)
  for (k in seq_along(groups)) {
    idx <- groups[[k]]
    S[idx, idx] <- 1
  }

  termino1 <- sum((3 * C - 2 * n) * S)
  termino2 <- sum(2 * (n - C))
  numerador <- termino1 + termino2
  denominador <- n * m*(2 * m - 1)
  QI <- numerador / denominador
  return(QI)
}

#' Robust Silhouette Calculation
#'
#' @param clusters Integer vector of cluster assignments
#' @param d_mat Distance matrix
#' @return Mean silhouette value
#' @export
calculate_silhouette_robust <- function(clusters, d_mat) {
  tryCatch({
    sil <- cluster::silhouette(clusters, dist = d_mat)
    return(mean(sil[, 3]))
  }, error = function(e) {
    n <- length(clusters)
    k <- length(unique(clusters))
    if (k < 2 || k >= n) return(NA)
    sil_vals <- numeric(n)
    for (i in 1:n) {
      cluster_i <- clusters[i]
      same_cluster <- which(clusters == cluster_i)
      if (length(same_cluster) > 1) {
        a_i <- mean(d_mat[i, same_cluster[same_cluster != i]])
      } else {
        a_i <- 0
      }
      b_i <- Inf
      for (cl in unique(clusters)) {
        if (cl != cluster_i) {
          other_cluster <- which(clusters == cl)
          mean_dist <- mean(d_mat[i, other_cluster])
          b_i <- min(b_i, mean_dist)
        }
      }
      sil_vals[i] <- (b_i - a_i) / max(a_i, b_i)
    }
    return(mean(sil_vals, na.rm = TRUE))
  })
}

#' Show partition with labels
#'
#' @param groups List of cluster memberships
#' @param labels Character vector of observation labels
#' @param method_name Name of the method
#' @return Invisible NULL, prints partition
#' @keywords internal
show_partition <- function(groups, labels = NULL, method_name = "Method") {
  cat(sprintf("\n=== Partition from %s ===\n", method_name))
  if (is.null(labels)) labels <- as.character(seq_along(groups[[1]]))
  for (i in seq_along(groups)) {
    idx <- groups[[i]]
    etiqueta <- if (!is.null(labels)) labels[idx] else idx
    cat(sprintf("  Cluster %d (size=%d): %s\n", i, length(idx),
                paste(etiqueta, collapse = ", ")))
  }
  invisible(NULL)
}

#' Internal: compute Gower distance matrix
#' @keywords internal
.compute_gower <- function(data) {
  as.matrix(cluster::daisy(data, metric = "gower"))
}