#' DBSCAN Clustering
#'
#' @param data A data.frame with categorical variables
#' @param eps_range Range of eps values to search
#' @param minPts_range Range of minPts values to search
#' @return List with results
#' @export
dbscan_clustering <- function(data, eps_range = seq(0.15, 0.55, by = 0.05),
                              minPts_range = 2:5) {
  d_gower <- .compute_gower(data)
  best_qi <- -Inf
  grupos_best <- NULL
  k_best <- NA
  n_noise_best <- NA

  for (eps in eps_range) {
    for (minPts in minPts_range) {
      res <- dbscan::dbscan(d_gower, eps = eps, minPts = minPts)
      n_clusters <- max(res$cluster)
      n_noise <- sum(res$cluster == 0)
      if (n_clusters >= 2) {
        grupos <- split(which(res$cluster > 0), res$cluster[res$cluster > 0])
        if (length(grupos) >= 2) {
          qi <- calculate_qi(data, grupos)
          coverage <- (nrow(data) - n_noise) / nrow(data)
          qi_adjusted <- qi * coverage
          if (qi_adjusted > best_qi) {
            best_qi <- qi_adjusted
            grupos_best <- grupos
            k_best <- n_clusters
            n_noise_best <- n_noise
          }
        }
      }
    }
  }

  list(
    method = "DBSCAN",
    k_opt = k_best,
    QI = best_qi,
    groups = grupos_best,
    noise_count = n_noise_best
  )
}

#' HDBSCAN Clustering
#'
#' @param data A data.frame with categorical variables
#' @param minPts_range Range of minPts values to search
#' @param seed Random seed for reproducibility
#' @return List with results
#' @export
hdbscan_clustering <- function(data, minPts_range = 2:10, seed = 42) {
  d_gower <- .compute_gower(data)
  best_qi <- -Inf
  grupos_best <- NULL
  k_best <- NA
  n_noise_best <- NA

  for (minPts in minPts_range) {
    set.seed(seed)
    res <- tryCatch({
      dbscan::hdbscan(d_gower, minPts = minPts)
    }, error = function(e) NULL)

    if (!is.null(res)) {
      n_clusters <- max(res$cluster)
      n_noise <- sum(res$cluster == 0)
      if (n_clusters >= 2) {
        grupos <- split(which(res$cluster > 0), res$cluster[res$cluster > 0])
        if (length(grupos) >= 2) {
          qi <- calculate_qi(data, grupos)
          coverage <- (nrow(data) - n_noise) / nrow(data)
          qi_adjusted <- qi * coverage
          if (qi_adjusted > best_qi) {
            best_qi <- qi_adjusted
            grupos_best <- grupos
            k_best <- n_clusters
            n_noise_best <- n_noise
          }
        }
      }
    }
  }

  list(
    method = "HDBSCAN",
    k_opt = k_best,
    QI = best_qi,
    groups = grupos_best,
    noise_count = n_noise_best
  )
}