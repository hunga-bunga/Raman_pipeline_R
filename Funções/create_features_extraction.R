#Extração de FEATURES 
create_features_dataset <- function(d_spectra_names, d_spectra_data, vector, smooth, folder1, folder2, ...) {
  
  # 1. Definir os nomes das colunas de características para cada valor do vetor
  peak_columns <- c('height', 'prominence', 'width', 'width_height', 'left', 'right', 'AUC')
  
  # Criar dinamicamente os nomes combinados (ex: "520_height", "520_prominence"...)
  new_columns <- c()
  for (val in sort(vector)) {
    for (col in peak_columns) {
      new_columns <- c(new_columns, paste0(val, "_", col))
    }
  }
  all_columns <- c(new_columns, "Label")
  
  # 2. Inicializar uma matriz vazia para armazenar os dados (linhas = amostras, colunas = características)
  num_samples <- length(d_spectra_names)
  num_features <- length(new_columns)
  data_matrix <- matrix(0, nrow = num_samples, ncol = num_features)
  colnames(data_matrix) <- new_columns
  labels <- character(num_samples)
  
  # 3. Loop por cada ficheiro/espectro (equivalente ao for file, raman_spectrum em...)
  for (idx in 1:num_samples) {
    file_name <- d_spectra_names[idx]
    raman_spectrum <- d_spectra_data[[idx]] # Assumindo lista de espectros
    
    # --- SIMULAÇÃO DO TEU PIPELINE INTERNO ---
    # Nota: Substitui estas duas linhas pelas chamadas reais do teu pacote/pipeline se necessário.
    # Aqui assume-se que extrais uma lista com os vetores de picos detetados.
    
    # Exemplo do que get_features devolve no Python:
    # properties <- list(Peaks = c(...), peak_heights = c(...), AUC = c(...))
    properties <- extract_spectrum_features(raman_spectrum, smooth, folder1, folder2, ...) 
    
    # Extrair a Label simulando o file.split('$')[0] do Python
    labels[idx] <- strsplit(file_name, "\\$")[[1]][1]
    
    # --- LOGICA ORIGINAL DE BUSCA DE PICOS (PEAK FEATURES) ---
    for (v in vector) {
      min_v <- v - 5
      max_v <- v + 5
      closest_peak_distance <- Inf
      closest_peak_index <- NULL
      
      # Procura o pico mais próximo dentro da janela [v-5, v+5]
      if (length(properties$Peaks) > 0) {
        for (i in 1:length(properties$Peaks)) {
          value <- properties$Peaks[i]
          
          if (value >= min_v && value <= max_v) {
            peak_distance <- abs(value - v)
            
            if (peak_distance < closest_peak_distance) {
              closest_peak_distance <- peak_distance
              closest_peak_index <- i
            }
          }
        }
      }
      
      # Se encontrou o pico, popula as colunas correspondentes na matriz
      if (!is.null(closest_peak_index)) {
        data_matrix[idx, paste0(v, "_height")]       <- properties$peak_heights[closest_peak_index]
        data_matrix[idx, paste0(v, "_prominence")]   <- properties$prominences[closest_peak_index]
        data_matrix[idx, paste0(v, "_width")]        <- properties$widths[closest_peak_index]
        data_matrix[idx, paste0(v, "_width_height")] <- properties$width_heights[closest_peak_index]
        data_matrix[idx, paste0(v, "_left")]         <- properties$left_minima[closest_peak_index]
        data_matrix[idx, paste0(v, "_right")]        <- properties$right_minima[closest_peak_index]
        data_matrix[idx, paste0(v, "_AUC")]          <- properties$AUC[closest_peak_index]
      } else {
        # Se não encontrou, os valores continuam a zero (já inicializados como 0 na matriz)
      }
    }
  }
  
  # 4. Juntar a matriz de dados com as Labels num data.frame estruturado
  df_features <- as.data.frame(data_matrix)
  df_features$Label <- as.factor(labels)
  
  return(df_features)
}
