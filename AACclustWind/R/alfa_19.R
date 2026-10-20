#' Internal AAC19 function (pure R implementation)
#'
#' @param data A data.frame or matrix with categorical variables
#' @param alpha Number of clusters desired
#' @return A list of integer vectors (cluster memberships)
#' @keywords internal
alfa19 <- function(data, alpha) {
  data_mat <- as.matrix(data)
  m <- nrow(data_mat)
  n <- ncol(data_mat)
  
  # VALIDATION: Dataset less than 3
  if (m < 3) {
    warning("Dataset very small. Trivial partition.")
    if (alpha == 1) return(list(1:m))
    if (alpha >= m) return(as.list(1:m))
    return(as.list(1:m))
  }
  
  # Funcion auxiliar SEGURA para fusionar clusters
  fusionar_clusters_seguro <- function(grupos, ab, B, dup) {
    if (length(dup) < 2) {
      return(list(grupos = grupos, ab = ab, A = length(grupos) + length(B[-ab])))
    }
    ak <- which(sapply(grupos, function(e) is.element(dup[1], e)))
    if (length(ak) < 2) {
      return(list(grupos = grupos, ab = ab, A = length(grupos) + length(B[-ab])))
    }
    if (length(dup) == 4 && length(ak) >= 2) {
      grupos[[ak[1]]] <- unique(unlist(grupos[ak]))
      grupos[[ak[2]]] <- NULL
      ab <- unlist(grupos)
      dup <- ab[!(!duplicated(ab) & rev(!duplicated(rev(ab))))]
      if (length(dup) >= 2) {
        ak <- which(sapply(grupos, function(e) is.element(dup[1], e)))
        if (length(ak) >= 2) {
          grupos[[ak[1]]] <- unique(unlist(grupos[ak]))
          grupos[[ak[2]]] <- NULL
        }
      }
    } else if (length(ak) >= 2) {
      az <- which(sapply(grupos, function(e) is.element(dup[1], e)))
      if (length(az) >= 2) {
        grupos[[az[1]]] <- unique(unlist(grupos[az]))
        grupos[[az[length(az)]]] <- NULL
      }
    }
    ab <- unlist(grupos)
    A <- length(grupos) + length(B[-ab])
    return(list(grupos = grupos, ab = ab, A = A))
  }
  
  # Matriz de coincidencia C
  C <- matrix(0, m, m)
  for (i in 1:m) {
    C[i, ] <- colSums(t(data_mat) == data_mat[i, ])
  }
  
  enlace <- C
  enlace <- enlace - diag(diag(enlace))
  grupos <- list()
  A <- m
  B <- 1:m
  suma <- 1
  
  # PRIMER WHILE: Construccion de clusters
  while (A > alpha && suma != 0) {
    k <- numeric(m)
    for (j in 1:m) {
      k[j] <- max(enlace[, j])
    }
    if (A == m) {
      ini <- which.max(k)
    } else {
      k[ab] <- 0
      if (sum(k) == 0) break
      ini <- which.max(k)
    }
    posicion <- which.max(enlace[ini, ])
    a <- enlace[ini, posicion]
    di <- which(enlace[ini, ] == a)
    enlace[ini, posicion] <- 0
    posicion_columna <- which.max(enlace[, posicion])
    b <- enlace[posicion_columna, posicion]
    enlace[ini, posicion] <- a
    
    if (a > b || a == b || a < b) {
      if (length(di) > 1) {
        if (m == A) {
          grupos[[1]] <- c(ini, posicion)
          ab <- unlist(grupos)
          largo <- length(B[-ab])
          A <- length(grupos) + largo
          enlace[, posicion] <- 0
          enlace[posicion, ] <- 0
        } else {
          grupos[[length(grupos) + 1]] <- c(ini, posicion)
          ab <- unlist(grupos)
          dup <- ab[!(!duplicated(ab) & rev(!duplicated(rev(ab))))]
          resultado <- fusionar_clusters_seguro(grupos, ab, B, dup)
          grupos <- resultado$grupos
          ab <- resultado$ab
          A <- resultado$A
          enlace[, posicion] <- 0
          enlace[posicion, ] <- 0
        }
        di2 <- which(enlace[ini, ] == a)
        for (j in 1:length(di2)) {
          if (A > alpha) {
            posicion <- which.max(enlace[ini, ])
            a2 <- enlace[ini, posicion]
            enlace[ini, posicion] <- 0
            posicion_columna <- which.max(enlace[, posicion])
            b2 <- enlace[posicion_columna, posicion]
            enlace[ini, posicion] <- a2
            if (a2 >= b2) {
              enlace[, posicion] <- 0
              enlace[posicion, ] <- 0
              grupos[[length(grupos) + 1]] <- c(ini, posicion)
              ab <- unlist(grupos)
              dup <- ab[!(!duplicated(ab) & rev(!duplicated(rev(ab))))]
              resultado <- fusionar_clusters_seguro(grupos, ab, B, dup)
              grupos <- resultado$grupos
              ab <- resultado$ab
              A <- resultado$A
            }
          }
        }
      } else {
        if (m == A) {
          grupos[[1]] <- c(ini, posicion)
          ab <- unlist(grupos)
          largo <- length(B[-ab])
          A <- length(grupos) + largo
          enlace[, posicion] <- 0
          enlace[posicion, ] <- 0
        } else {
          grupos[[length(grupos) + 1]] <- c(ini, posicion)
          enlace[, posicion] <- 0
          enlace[posicion, ] <- 0
          ab <- unlist(grupos)
          dup <- ab[!(!duplicated(ab) & rev(!duplicated(rev(ab))))]
          resultado <- fusionar_clusters_seguro(grupos, ab, B, dup)
          grupos <- resultado$grupos
          ab <- resultado$ab
          A <- resultado$A
        }
      }
    }
    
    k <- numeric(m)
    for (j in 1:m) {
      k[j] <- max(enlace[, j])
    }
    ab <- unlist(grupos)
    k[ab] <- 0
    largo <- length(B[-ab])
    if (sum(k) == 0 || largo == 0) {
      suma <- 0
    } else {
      suma <- sum(k)
    }
    k2 <- length(k[k != 0])
    if (k2 == 1) {
      suma <- 0
    }
  }
  
  # SEGUNDO WHILE: Optimizacion de clusters
  while (A > alpha && suma == 0) {
    ab3 <- length(unlist(grupos))
    if (ab3 == m) {
      af2 <- 1:length(grupos)
      combinaciones <- list()
      combi <- t(combn(length(af2), 2))
      for (an in 1:dim(combi)[1]) {
        combinaciones[[an]] <- combi[an, ]
      }
      Calidad_c <- rep(0, length(combinaciones))
      for (kp in 1:length(Calidad_c)) {
        grupos2 <- grupos
        grupos2[[length(grupos2) + 1]] <- unlist(grupos2[combinaciones[[kp]]])
        grupos2[combinaciones[[kp]]] <- NULL
        class_C <- numeric()
        for (nn in 1:length(grupos2)) {
          class_C[grupos2[[nn]]] <- nn
        }
        S <- outer(class_C, class_C, "==") * 1
        termino1 <- sum((3 * C - 2 * n) * S)
        termino2 <- sum(2 * (n - C))
        numerador <- termino1 + termino2
        denominador <- m^2 * n
        Calidad_c[kp] <- numerador / denominador
      }
      max_calidad <- which.max(Calidad_c)
      posi_calidad <- unlist(combinaciones[max_calidad])
      grupos[[posi_calidad[1]]] <- unlist(grupos[posi_calidad])
      grupos[posi_calidad[2]] <- NULL
      ab <- unlist(grupos)
      largo <- length(B[-ab])
      A <- length(grupos) + largo
    } else {
      ab <- unlist(grupos)
      largo2 <- B[-ab]
      grupos[[length(grupos) + 1]] <- largo2
      af2 <- 1:length(grupos)
      combinaciones <- list()
      combi <- t(combn(length(af2), 2))
      for (an in 1:dim(combi)[1]) {
        combinaciones[[an]] <- combi[an, ]
      }
      Calidad_c <- rep(0, length(combinaciones))
      for (kp in 1:length(Calidad_c)) {
        grupos2 <- grupos
        grupos2[[length(grupos2) + 1]] <- unlist(grupos2[combinaciones[[kp]]])
        grupos2[combinaciones[[kp]]] <- NULL
        class_C <- numeric()
        for (nn in 1:length(grupos2)) {
          class_C[grupos2[[nn]]] <- nn
        }
        S <- matrix(0, ncol = m, nrow = m)
        for (i in 1:m) {
          for (j in 1:m) {
            S[i, j] <- ifelse(class_C[i] == class_C[j], 1, 0)
          }
        }
        termino1 <- sum((3 * C - 2 * n) * S)
        termino2 <- sum(2 * (n - C))
        numerador <- termino1 + termino2
        denominador <- m^2 * n
        Calidad_c[kp] <- numerador / denominador
      }
      max_calidad <- which.max(Calidad_c)
      posi_calidad <- unlist(combinaciones[max_calidad])
      grupos[[posi_calidad[1]]] <- unlist(grupos[posi_calidad])
      grupos[posi_calidad[2]] <- NULL
      ab <- unlist(grupos)
      largo <- length(B[-ab])
      A <- length(grupos) + largo
    }
  }
  
  # ASIGNACION FINAL
  z <- list()
  z <- grupos
  ab <- unlist(grupos)
  ac <- B[-ab]
  if (length(grupos) == alpha) {
    z <- grupos
  } else {
    for (i in 1:length(ac)) {
      z[[length(grupos) + i]] <- ac[i]
    }
  }
  return(z)
}