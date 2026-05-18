# --- 1. FUNÇÃO AIRPLS (Sem funções auxiliares) ---
AIRPLS <- function(wavenumber, intensity, smoothing) {
  N <- length(intensity)
  
  # --- Construção da Matriz de Diferença de Ordem 2 (de raiz) ---
  D <- matrix(0, nrow = N - 2, ncol = N)
  for (i in 1:(N - 2)) {
    D[i, i]     <- 1
    D[i, i + 1] <- -2
    D[i, i + 2] <- 1
  }
  
  # Matriz de penalidade: P = lambda * D^T * D
  P <- smoothing * (t(D) %*% D)
  
  w <- rep(1, N)
  baseline_y <- rep(0, N)
  max_iter <- 500
  tol <- 0.001
  
  # --- Loop do Algoritmo airPLS ---
  for (i in 1:max_iter) {
    # Resolve o sistema linear (diag(w) + P) * baseline = w * intensity
    matriz_sistema <- diag(w) + P
    vetor_independente <- w * intensity
    baseline_y <- solve(matriz_sistema, vetor_independente)
    
    d <- intensity - baseline_y
    d_neg <- d[d < 0]
    
    if (length(d_neg) == 0) break
    
    sum_neg <- sum(abs(d_neg))
    if (i > 1 && (sum_neg / sum(abs(intensity))) < tol) break
    
    # Atualização dos pesos do airPLS
    w <- rep(0, N)
    w[d >= 0] <- 0
    w[d < 0] <- exp(i * d[d < 0] / sum_neg)
  }
  
  # Subtrai a linha de base encontrada
  baseline_subtracted_y <- intensity - baseline_y
  
  return(list(
    baseline = list(x = wavenumber, y = baseline_y),
    corrected = list(x = wavenumber, y = baseline_subtracted_y)
  ))
}
