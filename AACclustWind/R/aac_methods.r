#' AAC Sequential (Pure R implementation)
#'
#' @param data A data.frame or matrix with categorical variables
#' @param k_range Range of k values to evaluate (default: 1:nrow(data))
#' @param verbose Logical, print progress
#' @return List with k_opt, QI, groups, and all results
#' @export
aac_sequential <- function(data, k_range = NULL, verbose = FALSE) {
  data_mat <- as.matrix(data)
  m <- nrow(data_mat)
  if (is.null(k_range)) k_range <- 1:m
  
  qi_vals <- numeric(length(k_range))
  grupos_vals <- vector("list", length(k_range))
  
  for (i in seq_along(k_range)) {
    k <- k_range[i]
    if (k == 1) {
      grupos <- list(1:m)
    } else if (k >= m) {
      grupos <- as.list(1:m)
    } else {
      grupos <- alfa19(data_mat, k)
    }
    qi_vals[i] <- calculate_qi(data, grupos)
    grupos_vals[[i]] <- grupos
    if (verbose) cat(sprintf("  k=%d -> QI=%.4f\n", k, qi_vals[i]))
  }
  
  k_opt <- k_range[which.max(qi_vals)]
  list(
    method = "AAC_Sequential",
    k_opt = k_opt,
    QI = max(qi_vals),
    groups = grupos_vals[[which.max(qi_vals)]],
    qi_curve = stats::setNames(qi_vals, k_range),
    all_groups = grupos_vals
  )
}

#' AAC Parallel (R implementation with parLapply for Windows)
#'
#' @param data A data.frame or matrix
#' @param k_range Range of k values
#' @param n_cores Number of cores (default: detectCores() - 1)
#' @param verbose Logical
#' @return List with results
#' @export
aac_parallel <- function(data, k_range = NULL, n_cores = NULL, verbose = FALSE) {
  data_mat <- as.matrix(data)
  m <- nrow(data_mat)
  if (is.null(k_range)) k_range <- 1:m
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(cl, {
    library(AACclustWind)
  })
  parallel::clusterExport(cl, c("data_mat", "m"), envir = environment())
  
  resultados_par <- parallel::parLapply(cl, k_range, function(k) {
    if (k == 1) {
      grupos <- list(1:m)
    } else if (k >= m) {
      grupos <- as.list(1:m)
    } else {
      grupos <- alfa19(data_mat, k)
    }
    qi <- AACclustWind::calculate_qi(data_mat, grupos)
    list(k = k, QI = qi, grupos = grupos)
  })
  
  parallel::stopCluster(cl)
  
  qi_vals <- sapply(resultados_par, function(x) x$QI)
  k_opt <- k_range[which.max(qi_vals)]
  
  list(
    method = "AAC_Parallel",
    k_opt = k_opt,
    QI = max(qi_vals),
    groups = resultados_par[[which.max(qi_vals)]]$grupos,
    qi_curve = stats::setNames(qi_vals, k_range),
    n_cores = n_cores
  )
}

#' AAC Rcpp+Parallel (C++ implementation with parallel)
#'
#' @param data A data.frame or matrix
#' @param k_range Range of k values
#' @param n_cores Number of cores
#' @param verbose Logical
#' @return List with results
#' @export
aac_rcpp_parallel <- function(data, k_range = NULL, n_cores = NULL, verbose = FALSE) {
  data_mat <- as.matrix(data)
  m <- nrow(data_mat)
  n <- ncol(data_mat)
  if (is.null(k_range)) k_range <- 1:m
  if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  
  C_global <- calcular_matriz_C_cpp(data_mat)
  
  run_one <- function(k) {
    if (k == 1) {
      grupos <- list(1:m)
    } else if (k >= m) {
      grupos <- as.list(1:m)
    } else {
      grupos <- alfa19_cpp(data_mat, k)
    }
    qi <- calcular_QI_particion_cpp(C_global, grupos, m, n)
    list(k = k, QI = qi, grupos = grupos)
  }
  
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, {
      library(AACclustWind)
    })
    parallel::clusterExport(cl, c("data_mat", "m", "n", "C_global", "run_one"),
                            envir = environment())
    resultados_rcpp <- parallel::parLapply(cl, k_range, function(k) run_one(k))
    parallel::stopCluster(cl)
  } else {
    resultados_rcpp <- parallel::mclapply(k_range, run_one, mc.cores = n_cores)
  }
  
  qi_vals <- sapply(resultados_rcpp, function(x) x$QI)
  k_opt <- k_range[which.max(qi_vals)]
  
  list(
    method = "AAC_Rcpp_Parallel",
    k_opt = k_opt,
    QI = max(qi_vals),
    groups = resultados_rcpp[[which.max(qi_vals)]]$grupos,
    qi_curve = stats::setNames(qi_vals, k_range),
    n_cores = n_cores
  )
}