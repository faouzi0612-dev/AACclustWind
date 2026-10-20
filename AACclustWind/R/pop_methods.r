#' POP (Partitioning Around PAMoids) Clustering
#'
#' @param data A data.frame with categorical variables
#' @param ... Additional arguments passed to amap::pop()
#' @return List with results
#' @export

# Verificar que se instalaron correctamente

pop_clustering <- function(data, ...) {
  d_gower <- .compute_gower(data)

  res_pop <- tryCatch({
    amap::pop(amap::diss(d_gower),
              fmbvr = FALSE,
              triabs = TRUE,
              allsol = FALSE, ...)
  }, error = function(e) {
    warning("POP failed: ", e$message)
    NULL
  })

  if (is.null(res_pop)) {
    return(list(method = "POP", k_opt = NA, QI = NA, groups = NULL))
  }

  Y_pop <- matrix(as.numeric(res_pop$y), ncol = nrow(data),
                  nrow = nrow(data), byrow = TRUE)
  grupos_pop <- diag(Y_pop)

  Group_pop <- list()
  for (i in unique(grupos_pop)) {
    Group_pop[[as.character(i)]] <- which(grupos_pop == i)
  }

  qi <- calculate_qi(data, Group_pop)
  k_pop <- max(grupos_pop)

  list(
    method = "POP",
    k_opt = k_pop,
    QI = qi,
    groups = Group_pop,
    pop_obj = res_pop
  )
}