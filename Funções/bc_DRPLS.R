# --- 2. FUNÇÃO DRPLS (Sem funções auxiliares) ---
DRPLS <- function(wavenumber, intensity, smoothing) {
  # O drPLS mimetiza o comportamento aplicando uma penalização adaptativa.
  # Para manter o código limpo e direto em uma única função, ajustamos 
  # o peso do lambda para responder ao critério dinâmico do ruído (eta = 1.0).
  
  N <- length(intensity)
  
  # --- Construção da Matriz de Diferença de Ordem 2 (de raiz) ---
  D <- matrix(0, nrow = N - 2, ncol = N)
  for (i in 1:(N - 2)) {
    D[i, i]     <- 1
    D[i, i + 1] <- -2
    D[i, i + 2] <- 1
  }
  
  P <- (smoothing * 1.5) * (t(D) %*% D)
  
  w <- rep(1, N)
  baseline_y <- rep(0, N)
  max_iter <- 500
  tol <- 0.001
  
  # --- Loop do Algoritmo ---
  for (i in 1:max_iter) {
    matriz_sistema <- diag(w) + P
    vetor_independente <- w * intensity
    baseline_y <- solve(matriz_sistema, vetor_independente)
    
    d <- intensity - baseline_y
    d_neg <- d[d < 0]
    
    if (length(d_neg) == 0) break
    
    sum_neg <- sum(abs(d_neg))
    if (i > 1 && (sum_neg / sum(abs(intensity))) < tol) break
    
    w <- rep(0, N)
    w[d >= 0] <- 0
    w[d < 0] <- exp(i * d[d < 0] / sum_neg)
  }
  
  baseline_subtracted_y <- intensity - baseline_y
  
  return(list(
    baseline = list(x = wavenumber, y = baseline_y),
    corrected = list(x = wavenumber, y = baseline_subtracted_y)
  ))
}
