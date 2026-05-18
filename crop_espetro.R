# --- FUNÇÃO DE CROPPING ---
crop_espetro <- function(dataset, min_wn, max_wn) {
  # 1. Identificar índices (usa wavenumbers que agora existem no df)
  indices_manter <- which(dataset$wavenumbers >= min_wn & dataset$wavenumbers <= max_wn)
  
  if (length(indices_manter) == 0) {
    limites <- range(dataset$wavenumbers)
    stop(paste("Erro: Intervalo fora dos limites. O ficheiro vai de", 
               round(limites[1],1), "a", round(limites[2],1)))
  }
  
  # 2. Cortar
  dataset$data <- dataset$data[, indices_manter, drop = FALSE]
  dataset$wavenumbers <- dataset$wavenumbers[indices_manter]
  
  return(dataset)
}
