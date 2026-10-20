#' Hierarchical Clustering (Ward + Gower)
#'
#' @param data A data.frame with categorical variables
#' @param k_range Range of k to evaluate
#' @param method Linkage method (default: "ward.D2")
#' @return List with results
#' @export
hierarchical_clustering <- function(data, k_range = NULL, method = "ward.D2") {
  m <- nrow(data)
  if (is.null(k_range)) k_range <- 2:m

  d_gower <- .compute_gower(data)
  hc <- stats::hclust(stats::as.dist(d_gower), method = method)

  sil_vals <- c()
  k_vals <- c()

  for (k in k_range) {
    clusters <- stats::cutree(hc, k = k)
    if (length(unique(clusters)) == k) {
      sil_mean <- calculate_silhouette_robust(clusters, d_gower)
      if (!is.na(sil_mean)) {
        sil_vals <- c(sil_vals, sil_mean)
        k_vals <- c(k_vals, k)
      }
    }
  }

  if (length(sil_vals) > 0) {
    k_opt <- k_vals[which.max(sil_vals)]
    grupos <- split(1:m, stats::cutree(hc, k = k_opt))
  } else {
    k_opt <- min(5, m)
    grupos <- split(1:m, stats::cutree(hc, k = k_opt))
  }

  qi <- calculate_qi(data, grupos)

  list(
    method = "Hierarchical_Ward",
    k_opt = k_opt,
    QI = qi,
    groups = grupos,
    silhouette = max(sil_vals),
    hclust_obj = hc
  )
}

#' DIANA (Divisive Analysis) Clustering
#'
#' @param data A data.frame with categorical variables
#' @param k_range Range of k to evaluate
#' @return List with results
#' @export
diana_clustering <- function(data, k_range = NULL) {
  m <- nrow(data)
  if (is.null(k_range)) k_range <- 1:m

  d_gower <- .compute_gower(data)
  res_diana <- cluster::diana(d_gower)

  sil_vals <- c()
  k_vals <- c()

  for (k in k_range) {
    clusters <- stats::cutree(stats::as.hclust(res_diana), k = k)
    if (length(unique(clusters)) == k) {
      sil_mean <- calculate_silhouette_robust(clusters, d_gower)
      if (!is.na(sil_mean)) {
        sil_vals <- c(sil_vals, sil_mean)
        k_vals <- c(k_vals, k)
      }
    }
  }

  if (length(sil_vals) > 0) {
    k_opt <- k_vals[which.max(sil_vals)]
    grupos <- split(1:m, stats::cutree(stats::as.hclust(res_diana), k = k_opt))
  } else {
    k_opt <- min(5, m)
    grupos <- split(1:m, stats::cutree(stats::as.hclust(res_diana), k = k_opt))
  }

  qi <- calculate_qi(data, grupos)

  list(
    method = "DIANA",
    k_opt = k_opt,
    QI = qi,
    groups = grupos,
    silhouette = max(sil_vals),
    diana_obj = res_diana
  )
}