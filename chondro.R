setwd("C:/Users/rafae/Documents/Projeto")
library(specmine)
library(caret)


# --- FUNÇÃO PARA CARREGAR .npy ---
dataset_npy <- function(data_dir = getwd(), x_file, y_file, w_file) {
  
  read_npy_core <- function(file_name) {
    file_path <- file.path(data_dir, file_name)
    if (!file.exists(file_path)) stop(paste("Ficheiro não encontrado:", file_name))
    
    con <- file(file_path, "rb")
    on.exit(close(con))
    
    magic <- readBin(con, what = "raw", n = 6)
    if (!all(magic == as.raw(c(0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59)))) stop("Formato inválido.")
    
    readBin(con, what = "raw", n = 2)
    header_len <- readBin(con, what = "integer", n = 1, size = 2, endian = "little")
    header_raw_bytes <- readBin(con, what = "raw", n = header_len)
    header_text <- rawToChar(header_raw_bytes)
    
    shape_match <- regexec("'shape'\\s*:\\s*\\(([^)]+)\\)", header_text)
    shape_content <- regmatches(header_text, shape_match)[[1]][2]
    dims <- as.numeric(strsplit(gsub("\\s", "", shape_content), ",")[[1]])
    
    is_f8 <- grepl("'descr'\\s*:\\s*'[<>|]f8'", header_text)
    total_elements <- prod(dims)
    raw_data <- readBin(con, what = "numeric", n = total_elements, 
                        size = if(is_f8) 8 else 4, endian = "little")
    
    if (length(dims) == 1) return(as.vector(raw_data))
    return(matrix(raw_data, nrow = dims[1], ncol = dims[2], byrow = TRUE))
  }
  
  message("A carregar ficheiros especificados...")
  
  
  return(list(
    data = read_npy_core(x_file),
    y_train = as.vector(read_npy_core(y_file)),
    wavenumbers = as.numeric(read_npy_core(w_file))
  ))
}

#Alinhamento

align_silicon_peak <- function(wavenumber, intensity, target_peak = 520, search_window = 20) {
  
  # 1. Define os limites da janela de busca ao redor do pico alvo
  lower_bound <- target_peak - search_window
  upper_bound <- target_peak + search_window
  
  # 2. Cria a máscara lógica para filtrar a região do Silício
  mask <- (wavenumber > lower_bound) & (wavenumber < upper_bound)
  
  # Se nenhum ponto for encontrado na janela, interrompe com um erro
  if (!any(mask)) {
    stop("Nenhum dado encontrado na janela de busca especificada.")
  }
  
  wavenumber_si <- wavenumber[mask]
  intensity_si <- intensity[mask]
  
  # 3. Encontra a posição (wavenumber) onde ocorre a maior intensidade
  # which.max() equivale ao np.argmax() do Python
  peak_x <- wavenumber_si[which.max(intensity_si)]
  
  # 4. Computa o desvio (shift)
  shift <- target_peak - peak_x
  
  # 5. Corrige o eixo X inteiro somando o desvio
  aligned_wavenumber <- wavenumber + shift
  
  # Retorna uma lista contendo o novo eixo X corrigido e o valor do shift calculado
  return(list(
    aligned_wavenumber = aligned_wavenumber,
    shift = shift
  ))
}



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



#Função Despiking

dataset_despike <- function(dataset, ma = 10, threshold = 7) {
  
  # 1. Definição da lógica interna
  clean_single_spectrum <- function(y) {
    n <- length(y)
    med_y <- median(y)
    mad_y <- median(abs(y - med_y))
    
    # Proteção contra divisão por zero
    if (mad_y == 0) return(y)
    
    # Cálculo do Modified Z-Score
    y_z <- 0.6745 * (y - med_y) / mad_y
    spikes <- abs(y_z) > threshold
    spike_idx <- which(spikes)
    
    y_out <- y
    
    # Substituição dos spikes pela média da vizinhança
    for (i in spike_idx) {
      w <- seq(max(1, i - ma), min(n, i + ma))
      we <- w[!spikes[w]] # Usar apenas pontos que não são spikes para a média
      if (length(we) > 0) {
        y_out[i] <- mean(y[we])
      }
    }
    return(y_out)
  }
  
  message("A executar Despiking no dataset...")
  
  # 2. Aplicação à matriz inteira
  # t(apply(...)) garante que a estrutura [amostras x wavenumbers] se mantém
  dataset$data <- t(apply(dataset$data, 1, clean_single_spectrum))
  
  # Registo da operação para histórico
  cat("Concluído: ", sum(dataset$data != dataset$data), " potenciais spikes tratados.\n")
  
  return(dataset)
}


#Savitzky-golay (Denoising)

savitzky_golay(dataset, p.order, window, deriv = 0)

#Baseline corretion
baseline_correction(dataset, method = "als")

# --- FUNÇÃO airPLS OPTIMIZADA (Vetorizada, Esparsa e Alinhada ao Zero) ---
airPLS_fast_dataset <- function(dataset, lambda = 1e5, max_iter = 50) {
  
  if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
  library(Matrix)
  
  message("A calcular a baseline via airPLS Esparso (Alta Performance)...")
  
  X <- dataset$data
  N_amostras <- nrow(X)
  N_wavenumbers <- ncol(X)
  
  I <- Diagonal(N_wavenumbers)
  D <- diff(I, differences = 2)
  P <- lambda * t(D) %*% D
  
  X_corrected <- matrix(0, nrow = N_amostras, ncol = N_wavenumbers)
  tol <- 0.001
  
  for (j in 1:N_amostras) {
    intensity <- X[j, ]
    w <- rep(1, N_wavenumbers)
    baseline_y <- rep(0, N_wavenumbers)
    
    for (i in 1:max_iter) {
      Z <- Diagonal(x = w) + P
      baseline_y <- as.numeric(solve(Z, w * intensity))
      
      d <- intensity - baseline_y
      d_neg <- d[d < 0]
      if (length(d_neg) == 0) break
      
      sum_neg <- sum(abs(d_neg))
      if (i > 1 && (sum_neg / sum(abs(intensity))) < tol) break
      
      w <- rep(0, N_wavenumbers)
      w[d < 0] <- exp(i * d[d < 0] / sum_neg)
    }
    
    # [CORREÇÃO]: Subtrai o mínimo individual para forçar o espetro a começar no 0
    sinal_puro <- intensity - baseline_y
    X_corrected[j, ] <- sinal_puro - min(sinal_puro)
  }
  
  dataset$data <- X_corrected
  return(dataset)
}

# --- 2. FUNÇÃO DRPLS (Sem funções auxiliares) ---
DRPLS <- function(wavenumber, intensity, smoothing) {

  
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

#Normalização

# --- FUNÇÃO DE NORMALIZAÇÃO EXPANDIDA ---
raman_norm <- function(dataset, method = "min-max", target_peak = 520, search_window = 20) {
  
  # Extrair a matriz de dados (amostras nas linhas, wavenumbers nas colunas)
  matriz_dados <- dataset$data
  wavenumber <- dataset$wavenumbers
  
  if (method == "min-max") {
    message("A aplicar Normalização Min-Max...")
    # (x - min) / (max - min) para cada linha
    norm_data <- t(apply(matriz_dados, 1, function(x) {
      (x - min(x)) / (max(x) - min(x))
    }))
    
  } else if (method == "L2") {
    message("A aplicar Normalização L2 (Vetorial)...")
    # x / sqrt(sum(x^2)) para cada linha
    norm_data <- t(apply(matriz_dados, 1, function(x) {
      norma_l2 <- sqrt(sum(x^2))
      if(norma_l2 == 0) return(x) # Evitar divisão por zero
      x / norma_l2
    }))
    
  } else if (method == "silicon") {
    message("A aplicar Normalização pelo Pico de Silício (Eixo Y)...")
    
    # 1. Define os limites da janela de busca ao redor do pico alvo
    lower_bound <- target_peak - search_window
    upper_bound <- target_peak + search_window
    
    # 2. Cria a máscara lógica usando os wavenumbers do dataset
    mask <- (wavenumber > lower_bound) & (wavenumber < upper_bound)
    
    if (!any(mask)) {
      stop("Erro: Nenhum dado encontrado na janela de busca especificada para o Silício.")
    }
    
    # 3. Processa cada amostra (linha) para normalizar pelo valor máximo do pico de Si
    norm_data <- t(apply(matriz_dados, 1, function(x) {
      # Isola as intensidades da janela daquela amostra específica
      intensity_si <- x[mask]
      # Encontra o valor máximo de intensidade no pico do Silício
      silica_peak_value <- max(intensity_si, na.rm = TRUE)
      
      if (silica_peak_value == 0) return(x) # Evitar divisão por zero
      
      # Divide o espectro inteiro pelo valor desse pico
      x / silica_peak_value
    }))
    
  } else {
    stop("Método inválido. Escolha entre 'min-max', 'L2' ou 'silicon'.")
  }
  
  # Atualiza o dataset com os novos dados e retornar
  dataset$data <- norm_data
  return(dataset)
}

#PCA E GRAFICO
raman_pca_plot <- function(dataset) {
  message("A calcular o PCA e a gerar o gráfico de clusters...")
  
  # 1. Executar o PCA na matriz de dados
  # scale. = TRUE é importante se os dados não estiverem perfeitamente normalizados
  pca_res <- prcomp(dataset$data, center = TRUE, scale. = FALSE)
  
  # 2. Calcula a variância explicada pelas duas primeiras componentes (PC1 e PC2)
  var_explicada <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  lbl_pc1 <- paste0("PC1 (", round(var_explicada[1], 1), "%)")
  lbl_pc2 <- paste0("PC2 (", round(var_explicada[2], 1), "%)")
  
  # 3. Cria um data frame com as coordenadas dos pontos e as etiquetas originais
  df_pca <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Classe = as.factor(dataset$y_train)
  )
  
  # 4. Desenhaa o gráfico de dispersão (Scatter Plot) nativo do R
  # Definir cores bonitas para as tuas 5 classes
  cores_classes <- c("#e63946", "#457b9d", "#1d3557", "#2a9d8f", "#f4a261")
  
  plot(df_pca$PC1, df_pca$PC2, 
       col = cores_classes[df_pca$Classe], 
       pch = 19, cex = 1.2,
       xlab = lbl_pc1, ylab = lbl_pc2,
       main = "Análise de Componentes Principais (PCA) - Clusters Raman")
  
  # Adicionar uma grelha de fundo
  grid(col = "gray90", lty = "dotted")
  
  # Adicionar a legenda com os nomes reais das tuas classes
  legend("topright", legend = levels(df_pca$Classe), 
         col = cores_classes, pch = 19, bty = "n", cex = 0.9)
  
  return(pca_res)
}


#SUPERVISED MACHINE LEARNING

ML_raman <- function(dataset, method = "rf", train_prop = 0.7) {
  df_ml <- as.data.frame(dataset$data)
  df_ml$class <- as.factor(dataset$y_train)
  
  set.seed(123)
  # Divisão por Paciente se disponível para evitar Data Leakage
  if (!is.null(dataset$patient_ids)) {
    message("A realizar a divisão estrita por PACIENTE no ML Clássico...")
    pacientes_unicos <- unique(dataset$patient_ids)
    pacientes_treino <- sample(pacientes_unicos, size = round(length(pacientes_unicos) * train_prop))
    trainIndex <- which(dataset$patient_ids %in% pacientes_treino)
    
    train_data <- df_ml[trainIndex, ]
    test_data  <- df_ml[-trainIndex, ]
    
    # Validação cruzada que também respeita o agrupamento de pacientes nos folds
    fitControl <- trainControl(method = "cv", number = 5,
                               index = groupKFold(dataset$patient_ids[trainIndex], k = 5))
  } else {
    warning("Aviso: 'patient_ids' não encontrado. A usar divisão por espetro (Risco de Data Leakage).")
    trainIndex <- createDataPartition(df_ml$class, p = train_prop, list = FALSE)
    
    train_data <- df_ml[trainIndex, ]
    test_data  <- df_ml[-trainIndex, ]
    
    fitControl <- trainControl(method = "cv", number = 5)
  }
  
  # --- [CORREÇÃO CRÍTICA: CÁLCULO DOS PESOS PARA O BALANCEAMENTO] ---
  # Calcula a frequência de cada classe no treino e inverte o peso
  frequencia_classes <- table(train_data$class)
  pesos_linhas <- 1 / as.numeric(frequencia_classes[train_data$class])
  
  message(paste("A treinar modelo:", method, "..."))
  
  # O caret utiliza o argumento 'weights' para penalizar erros nas classes minoritárias
  modelo <- switch(method,
                   "svm" = train(class ~ ., data = train_data, method = "svmRadial", 
                                 trControl = fitControl, weights = pesos_linhas),
                   "rf"  = train(class ~ ., data = train_data, method = "rf", 
                                 trControl = fitControl, weights = pesos_linhas),
                   "gbm" = train(class ~ ., data = train_data, method = "gbm", 
                                 trControl = fitControl, verbose = FALSE),
                   stop("Método inválido! Escolha 'svm', 'rf' ou 'gbm'.")
  )
  
  # Avaliação automática no conjunto de teste independente
  predicoes <- predict(modelo, newdata = test_data)
  cm <- confusionMatrix(predicoes, test_data$class)
  
  # Imprime os resultados na consola para conseguires analisar imediatamente
  cat("\n==================================================")
  cat(paste("\n[RESULTADOS DE VALIDACÃO -", toupper(method), "]:\n"))
  print(cm)
  cat("==================================================\n")
  
  return(modelo)
}

#Extração de FEATURES 
create_features_dataset <- function(d_spectra_names, properties_list, vector) {
  
  peak_columns <- c('height', 'prominence', 'width', 'width_height', 'left', 'right', 'AUC')
  
  # Criar dinamicamente os nomes das colunas combinadas
  new_columns <- c()
  for (val in sort(vector)) {
    for (col in peak_columns) {
      new_columns <- c(new_columns, paste0(val, "_", col))
    }
  }
  
  # Inicializar a matriz de dados e o vetor de labels
  num_samples <- length(d_spectra_names)
  num_features <- length(new_columns)
  data_matrix <- matrix(0, nrow = num_samples, ncol = num_features)
  colnames(data_matrix) <- new_columns
  labels <- character(num_samples)
  
  # Loop para processar cada amostra com os teus dados reais
  for (idx in 1:num_samples) {
    file_name <- d_spectra_names[idx]
    properties <- properties_list[[idx]] 
    
    # Extrai a Label a partir do nome do ficheiro (tudo antes do '$')
    labels[idx] <- strsplit(file_name, "\\$")[[1]][1]
    
    # Lógica de busca de picos dentro da janela de tolerância (+- 5 cm^-1)
    for (v in vector) {
      min_v <- v - 5
      max_v <- v + 5
      closest_peak_distance <- Inf
      closest_peak_index <- NULL
      
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
      
      # Se encontrar o pico real, preenche as respetivas características na matriz
      if (!is.null(closest_peak_index)) {
        data_matrix[idx, paste0(v, "_height")]       <- properties$peak_heights[closest_peak_index]
        data_matrix[idx, paste0(v, "_prominence")]   <- properties$prominences[closest_peak_index]
        data_matrix[idx, paste0(v, "_width")]        <- properties$widths[closest_peak_index]
        data_matrix[idx, paste0(v, "_width_height")] <- properties$width_heights[closest_peak_index]
        data_matrix[idx, paste0(v, "_left")]         <- properties$left_minima[closest_peak_index]
        data_matrix[idx, paste0(v, "_right")]        <- properties$right_minima[closest_peak_index]
        data_matrix[idx, paste0(v, "_AUC")]          <- properties$AUC[closest_peak_index]
      }
    }
  }
  
  # Converte para data.frame e define a Label como factor para o Machine Learning
  df_result <- as.data.frame(data_matrix)
  df_result$Label <- as.factor(labels)
  
  return(df_result)
}

ML_raman_cnn1d <- function(dataset, train_prop = 0.7, epochs = 15, batch_size = 32) {
  if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
  library(caret)
  if (!requireNamespace("keras3", quietly = TRUE)) install.packages("keras3")
  library(keras3)
  
  message("A preparar os dados para a CNN 1D...")
  X <- dataset$data
  y <- as.numeric(as.factor(dataset$y_train)) - 1
  num_classes <- length(unique(y))
  
  set.seed(123)
  # Divisão por Paciente se disponível, caso contrário por espetro
  if (!is.null(dataset$patient_ids)) {
    message("A realizar a divisão estrita por PACIENTE na CNN...")
    pacientes_unicos <- unique(dataset$patient_ids)
    pacientes_treino <- sample(pacientes_unicos, size = round(length(pacientes_unicos) * train_prop))
    trainIndex <- which(dataset$patient_ids %in% pacientes_treino)
  } else {
    warning("Aviso: 'patient_ids' não encontrado. A usar divisão por espetro (Risco de Data Leakage).")
    trainIndex <- createDataPartition(as.factor(y), p = train_prop, list = FALSE)
  }
  
  X_train <- X[trainIndex, , drop = FALSE]
  y_train <- y[trainIndex]
  X_test  <- X[-trainIndex, , drop = FALSE]
  y_test  <- y[-trainIndex]
  
  X_train_reshaped <- array(X_train, dim = c(nrow(X_train), ncol(X_train), 1))
  X_test_reshaped  <- array(X_test, dim = c(nrow(X_test), ncol(X_test), 1))
  
  y_train_categorical <- to_categorical(y_train, num_classes)
  y_test_categorical  <- to_categorical(y_test, num_classes)
  
  # --- CÁLCULO DE CLASS WEIGHTS ---
  tab_classes <- table(y_train)
  total_amostras <- sum(tab_classes)
  pesos_vetor <- total_amostras / (num_classes * as.numeric(tab_classes))
  class_weights_list <- as.list(pesos_vetor)
  names(class_weights_list) <- as.character(0:(num_classes - 1))
  
  input_shape <- c(ncol(X_train), 1)
  
  # Arquitetura Otimizada com Dropout Ajustado
  model <- keras_model_sequential() %>%
    layer_conv_1d(filters = 16, kernel_size = 5, activation = 'relu', 
                  kernel_regularizer = regularizer_l2(0.001), input_shape = input_shape) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_dropout(rate = 0.25) %>% 
    
    layer_conv_1d(filters = 32, kernel_size = 3, activation = 'relu',
                  kernel_regularizer = regularizer_l2(0.001)) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_dropout(rate = 0.3) %>% 
    
    layer_flatten() %>%
    layer_dense(units = 32, activation = 'relu', kernel_regularizer = regularizer_l2(0.001)) %>%
    layer_dropout(rate = 0.3) %>% 
    layer_dense(units = num_classes, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(learning_rate = 0.001),
    metrics = c('accuracy')
  )
  
  early_stop <- callback_early_stopping(
    monitor = "val_loss", 
    patience = 7, # Aumentado para dar mais tempo ao modelo balanceado para estabilizar
    restore_best_weights = TRUE
  )
  
  history <- model %>% fit(
    X_train_reshaped, y_train_categorical,
    epochs = epochs,
    batch_size = batch_size,
    validation_data = list(X_test_reshaped, y_test_categorical),
    callbacks = list(early_stop),
    class_weight = class_weights_list, # Injeção de pesos aplicada
    verbose = 1
  )
  
  y_pred_probs <- model %>% predict(X_test_reshaped)
  y_pred       <- apply(y_pred_probs, 1, which.max) - 1
  
  classes_originais <- levels(as.factor(dataset$y_train))
  y_test_factored   <- factor(classes_originais[y_test + 1], levels = classes_originais)
  y_pred_factored   <- factor(classes_originais[y_pred + 1], levels = classes_originais)
  
  cm <- confusionMatrix(y_pred_factored, y_test_factored)
  
  return(list(
    model = model,
    history = history,
    confusion_matrix = cm
  ))
}

#Wavelet e Fourier

raman_transform_wavelet <- function(dataset, level = 1) {
  message("A aplicar Wavelet Transform (Haar DWT)...")
  matriz_dados <- dataset$data
  
  # Função interna para aplicar DWT Haar a um vetor de raiz
  dwt_haar_1d <- function(x) {
    n <- length(x)
    # Se o comprimento for ímpar, ajustamos o último elemento
    if (n %% 2 != 0) x <- c(x, x[n])
    
    # Filtro Haar de raiz
    approx <- (x[seq(1, n, by=2)] + x[seq(2, n, by=2)]) / sqrt(2)
    detail <- (x[seq(1, n, by=2)] - x[seq(2, n, by=2)]) / sqrt(2)
    
    return(c(approx, detail))
  }
  
  # Aplica o nível de decomposição desejado
  transform_matrix <- t(apply(matriz_dados, 1, function(row) {
    res <- row
    for(l in 1:level) {
      # Decompõe apenas a parte de aproximação a cada nível
      sub_n <- length(res) / (2^(l-1))
      res[1:sub_n] <- dwt_haar_1d(res[1:sub_n])
    }
    return(res)
  }))
  
  dataset$data <- transform_matrix
  return(dataset)
}

raman_transform_fourier <- function(dataset) {
  message("A aplicar Fourier Transform (Power Spectrum)...")
  matriz_dados <- dataset$data
  
  power_spectrum <- t(apply(matriz_dados, 1, function(x) {
    # Aplica a Transformada Rápida de Fourier nativa do R
    fft_res <- fft(x)
    # Calcula a magnitude ao quadrado de cada componente de frequência (Power)
    # Retorna apenas a primeira metade devido à simetria do sinal real
    N <- length(x)
    mag <- (Mod(fft_res)^2) / N
    return(mag[1:floor(N/2)])
  }))
  
  dataset$data <- power_spectrum
  # Ajusta o vetor de frequências correspondente
  dataset$wavenumbers <- seq(0, 0.5, length.out = ncol(power_spectrum))
  
  return(dataset)
}

#Radar Gráfico/Visualização

raman_radar_plot <- function(dataset, sample_idx = 1) {
  message("A aplicar Transformação Radar e a gerar gráfico (Coordenadas Polares)...")
  
  matriz_dados <- dataset$data
  wavenumbers <- dataset$wavenumbers
  N_wavenumbers <- length(wavenumbers)
  
  # 1. Mapeia os wavenumbers uniformemente em ângulos de 0 a 2*pi
  theta <- seq(0, 2 * pi, length.out = N_wavenumbers)
  
  # 2. Converte as intensidades de todas as amostras para X e Y polares
  radar_xy_list <- list()
  for(idx in 1:nrow(matriz_dados)) {
    r <- matriz_dados[idx, ]
    
    # Normaliza o raio entre 0 e 1 para manter a proporção do círculo
    r_norm <- (r - min(r)) / (max(r) - min(r))
    if(any(is.na(r_norm))) r_norm <- rep(0, length(r))
    
    x_polar <- r_norm * cos(theta)
    y_polar <- r_norm * sin(theta)
    
    radar_xy_list[[idx]] <- data.frame(X_polar = x_polar, Y_polar = y_polar)
  }
  
  # 3. CONSTRUÇÃO AUTOMÁTICA DO GRÁFICO (Estilo da imagem)
  # Extrai os dados da amostra escolhida para o plot
  dados_plot <- radar_xy_list[[sample_idx]]
  
  # Configura a janela gráfica com proporção perfeitamente quadrada (asp = 1)
  plot(NULL, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), asp = 1,
       xlab = "", ylab = "", axes = FALSE,
       main = paste("Class Control - Average 3D Raman Spectra Projection\n(Amostra", sample_idx, ")"))
  
  # Desenha as linhas de grelha circulares (anéis concêntricos)
  valores_raio <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  for(raio in valores_raio) {
    # Desenha o círculo gerando pontos de 0 a 2*pi
    pts_circulo <- seq(0, 2*pi, length.out = 200)
    lines(raio * cos(pts_circulo), raio * sin(pts_circulo), col = "gray85", lty = 1)
  }
  
  # Desenha o círculo exterior principal (a borda preta fina)
  lines(1.05 * cos(pts_circulo), 1.05 * sin(pts_circulo), col = "black", lwd = 1.2)
  
  # Coloca os marcadores de Wavenumbers ao redor do círculo (como na imagem)
  # Seleciona 8 pontos cardeais/colaterais para legendar
  idx_labels <- round(seq(1, N_wavenumbers, length.out = 9))[1:8]
  for(i in idx_labels) {
    # Calcula a posição do texto ligeiramente fora do círculo (raio = 1.15)
    ang <- theta[i]
    txt_x <- 1.18 * cos(ang)
    txt_y <- 1.18 * sin(ang)
    
    # Arredonda o valor do wavenumber para mostrar na borda
    text(txt_x, txt_y, labels = round(wavenumbers[i]), cex = 0.8, col = "black")
  }
  
  # Preenche a área do espetro com a cor sólida de fundo (Azul Escuro da imagem)
  polygon(dados_plot$X_polar, dados_plot$Y_polar, col = "#0d1b2a", border = NA)
  
  # Desenha o contorno vibrante por cima (Efeito degradé/fogo da imagem)
  # Fazemos múltiplas camadas transparentes para dar o efeito de brilho (glow)
  polygon(dados_plot$X_polar, dados_plot$Y_polar, col = rgb(0.9, 0.1, 0.4, 0.4), border = NA)
  polygon(dados_plot$X_polar * 0.95, dados_plot$Y_polar * 0.95, col = rgb(1, 0.6, 0, 0.3), border = NA)
  
  # Desenha a linha de contorno final externa
  lines(dados_plot$X_polar, dados_plot$Y_polar, col = "#e63946", lwd = 1.5)
  
  # Adiciona um ponto azul brilhante no centro exato (origem)
  points(0, 0, col = "#00b4d8", pch = 16, cex = 1.2)
  
  # Salva o resultado no formato do pipeline e retorna
  dataset$data_polar <- radar_xy_list
  return(dataset)
}


raman_radar_class_averages <- function(dataset) {
  message("A calcular médias por classe e a gerar projeções polares (Radar Plot)...")
  
  matriz_dados <- dataset$data
  wavenumbers <- dataset$wavenumbers
  classes <- as.factor(dataset$y_train)
  unique_classes <- levels(classes)
  N_wavenumbers <- length(wavenumbers)
  
  # Ângulos circulares para o mapeamento completo do espetro (0 a 2*pi)
  theta <- seq(0, 2 * pi, length.out = N_wavenumbers)
  
  # Configurar o layout gráfico para mostrar todas as classes lado a lado numa linha
  # Exemplo: com 5 classes, divide a janela em 1 linha e 5 colunas
  par(mfrow = c(1, length(unique_classes)))
  
  # Cores vibrantes para os contornos de cada classe
  cores_linhas <- c("#e63946", "#457b9d", "#1d3557", "#2a9d8f", "#f4a261")
  
  for (i in 1:length(unique_classes)) {
    classe_atual <- unique_classes[i]
    
    # Isolar as amostras pertencentes a esta classe e calcular a média por coluna
    indices_classe <- which(classes == classe_atual)
    espetro_medio <- colMeans(matriz_dados[indices_classe, , drop = FALSE])
    
    # Normalizar o raio entre 0 e 1 baseado nos limites globais do espetro
    r_norm <- (espetro_medio - min(matriz_dados)) / (max(matriz_dados) - min(matriz_dados))
    
    x_polar <- r_norm * cos(theta)
    y_polar <- r_norm * sin(theta)
    
    # Desenhar o gráfico polar quadrado
    plot(NULL, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), asp = 1,
         xlab = "", ylab = "", axes = FALSE,
         main = paste("Classe:", classe_atual))
    
    # Grelha interna de anéis circulares
    for(raio in c(0.2, 0.4, 0.6, 0.8, 1.0)) {
      pts_circulo <- seq(0, 2*pi, length.out = 100)
      lines(raio * cos(pts_circulo), raio * sin(pts_circulo), col = "gray90", lty = 1)
    }
    lines(1.05 * cos(pts_circulo), 1.05 * sin(pts_circulo), col = "black", lwd = 1)
    
    # Polígono preenchido com fundo escuro e glow transparente
    polygon(x_polar, y_polar, col = "#0d1b2a", border = NA)
    polygon(x_polar, y_polar, col = rgb(0.9, 0.1, 0.4, 0.2), border = NA)
    
    # Linha de contorno com a cor específica daquela classe
    lines(x_polar, y_polar, col = cores_linhas[i], lwd = 2)
    
    # Centro
    points(0, 0, col = "#00b4d8", pch = 16, cex = 1)
  }
  
  # Resetar o layout de ecrã do RStudio para o padrão 1x1
  par(mfrow = c(1, 1))
}

visualizar_raw_vs_pre <- function(df_raw, df_processado, metodo_norm = "L2") {
  # Configura a janela gráfica para 2 gráficos verticais
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  
  # --- 1. GRÁFICO: DADOS ORIGINAIS (RAW) ---
  matriz_raw <- df_raw$data
  wavenumbers_raw <- df_raw$wavenumbers
  
  # Calcular métricas estatísticas dos dados brutos
  media_raw <- colMeans(matriz_raw, na.rm = TRUE)
  min_raw   <- apply(matriz_raw, 2, min, na.rm = TRUE)
  max_raw   <- apply(matriz_raw, 2, max, na.rm = TRUE)
  
  # Gráfico base vazio com os limites corretos do Raw
  plot(NULL, xlim = range(wavenumbers_raw), ylim = range(matriz_raw, na.rm = TRUE),
       xlab = "Wavenumber (cm-1)", ylab = "Intensidade (U.A.)",
       main = paste("Dataset ORIGINAL (Média Geral + Banda Total, N =", nrow(matriz_raw), ")"))
  grid(col = "gray85")
  
  # Desenhar a área de variabilidade (sombra cinzenta)
  polygon(c(wavenumbers_raw, rev(wavenumbers_raw)), c(max_raw, rev(min_raw)),
          col = rgb(0.8, 0.8, 0.8, alpha = 0.4), border = NA)
  
  # Desenhar a linha da média original a vermelho
  lines(wavenumbers_raw, media_raw, col = "red", lwd = 2.5)
  
  
  # --- 2. GRÁFICO: DADOS PRÉ-PROCESSADOS (DINÂMICO) ---
  matriz_proc <- df_processado$data
  wavenumbers_proc <- df_processado$wavenumbers
  
  # Calcular métricas estatísticas dos dados limpos
  media_proc <- colMeans(matriz_proc, na.rm = TRUE)
  min_proc   <- apply(matriz_proc, 2, min, na.rm = TRUE)
  max_proc   <- apply(matriz_proc, 2, max, na.rm = TRUE)
  
  # Definir cores e rótulos baseados no método escolhido
  if (toupper(metodo_norm) == "L2") {
    cor_linha  <- "darkblue"
    cor_sombra <- rgb(0.7, 0.8, 0.9, alpha = 0.5) # Sombra azulada para L2
    label_y    <- "Intensidade Normalizada L2"
    titulo_pos <- paste("Dataset PRÉ-PROCESSADO (Média Geral + Banda L2, N =", nrow(matriz_proc), ")")
  } else if (toupper(metodo_norm) == "MINMAX" || toupper(metodo_norm) == "MIN-MAX") {
    cor_linha  <- "darkgreen"
    cor_sombra <- rgb(0.7, 0.9, 0.7, alpha = 0.5) # Sombra esverdeada para MinMax
    label_y    <- "Intensidade Normalizada Min-Max [0, 1]"
    titulo_pos <- paste("Dataset PRÉ-PROCESSADO (Média Geral + Banda Min-Max, N =", nrow(matriz_proc), ")")
  } else {
    # Fallback genérico caso passes outro nome
    cor_linha  <- "black"
    cor_sombra <- rgb(0.6, 0.6, 0.6, alpha = 0.5)
    label_y    <- paste("Intensidade Normalizada", metodo_norm)
    titulo_pos <- paste("Dataset PRÉ-PROCESSADO (Média Geral + Banda", metodo_norm, ", N =", nrow(matriz_proc), ")")
  }
  
  # Desenhar o gráfico processado com as variáveis dinâmicas
  plot(NULL, xlim = range(wavenumbers_proc), ylim = range(matriz_proc, na.rm = TRUE),
       xlab = "Wavenumber (cm-1)", ylab = label_y, main = titulo_pos)
  grid(col = "gray85")
  
  # Desenhar a área de variabilidade pós-processamento
  polygon(c(wavenumbers_proc, rev(wavenumbers_proc)), c(max_proc, rev(min_proc)),
          col = cor_sombra, border = NA)
  
  # Desenhar a linha da média processada
  lines(wavenumbers_proc, media_proc, col = cor_linha, lwd = 2.5)
  
  # Restaurar definições de ecrã padrão do R
  par(mfrow = c(1, 1))
}

ML_raman_cnn2d_radar <- function(dataset, img_size = 64, train_prop = 0.7, epochs = 20, batch_size = 32) {
  if (!requireNamespace("keras3", quietly = TRUE)) install.packages("keras3")
  library(keras3)
  library(caret)
  
  message("A transformar espetros radar em imagens 2D matriciais...")
  
  matriz_dados <- dataset$data
  wavenumbers <- dataset$wavenumbers
  N_amostras <- nrow(matriz_dados)
  N_wavenumbers <- length(wavenumbers)
  theta <- seq(0, 2 * pi, length.out = N_wavenumbers)
  
  # Inicializar o array de imagens [Amostras, Altura, Largura, Canais]
  X_images <- array(0, dim = c(N_amostras, img_size, img_size, 1))
  
  # Gerar a matriz de imagem para cada amostra de forma eficiente
  for (idx in 1:N_amostras) {
    r <- matriz_dados[idx, ]
    # Normalizar o raio entre 0 e 1
    r_norm <- (r - min(r)) / (max(r) - min(r) + 1e-8)
    
    # Coordenadas polares mapeadas para o espaço [-1, 1]
    x_polar <- r_norm * cos(theta)
    y_polar <- r_norm * sin(theta)
    
    # Mapear coordenadas [-1, 1] para os índices da matriz [1, img_size]
    col_idx <- round(((x_polar + 1) / 2) * (img_size - 1)) + 1
    row_idx <- round(((y_polar + 1) / 2) * (img_size - 1)) + 1
    
    # Desenhar as linhas do radar preenchendo a matriz de pixéis
    for (i in 1:length(row_idx)) {
      X_images[idx, row_idx[i], col_idx[i], 1] <- 1.0 # Pixel ativo onde passa o sinal
    }
  }
  
  # Preparar as etiquetas (Labels)
  y_factors <- as.factor(dataset$y_train)
  classes_originais <- levels(y_factors)
  num_classes <- length(classes_originais)
  y <- as.numeric(y_factors) - 1
  
  # Divisão dos dados (Treino / Teste)
  set.seed(123)
  if (!is.null(dataset$patient_ids)) {
    pacientes_unicos <- unique(dataset$patient_ids)
    pacientes_treino <- sample(pacientes_unicos, size = round(length(pacientes_unicos) * train_prop))
    trainIndex <- which(dataset$patient_ids %in% pacientes_treino)
  } else {
    trainIndex <- createDataPartition(as.factor(y), p = train_prop, list = FALSE)
  }
  
  X_train <- X_images[trainIndex, , , , drop = FALSE]
  y_train <- y[trainIndex]
  X_test  <- X_images[-trainIndex, , , , drop = FALSE]
  y_test  <- y[-trainIndex]
  
  y_train_cat <- to_categorical(y_train, num_classes)
  y_test_cat  <- to_categorical(y_test, num_classes)
  
  # Arquitetura da CNN 2D para Reconhecimento de Imagem Radar
  model <- keras_model_sequential() %>%
    layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = 'relu', input_shape = c(img_size, img_size, 1)) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_dropout(rate = 0.25) %>%
    
    layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = 'relu') %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_dropout(rate = 0.25) %>%
    
    layer_flatten() %>%
    layer_dense(units = 64, activation = 'relu') %>%
    layer_dropout(rate = 0.4) %>%
    layer_dense(units = num_classes, activation = 'softmax')
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(learning_rate = 0.001),
    metrics = c('accuracy')
  )
  
  # Treino do modelo
  message("A treinar a CNN 2D com as imagens geradas...")
  history <- model %>% fit(
    X_train, y_train_cat,
    epochs = epochs,
    batch_size = batch_size,
    validation_data = list(X_test, y_test_cat),
    verbose = 1
  )
  
  # Predições e Matriz de Confusão
  y_pred_probs <- model %>% predict(X_test)
  y_pred       <- apply(y_pred_probs, 1, which.max) - 1
  
  y_test_factored <- factor(classes_originais[y_test + 1], levels = classes_originais)
  y_pred_factored <- factor(classes_originais[y_pred + 1], levels = classes_originais)
  
  cm <- confusionMatrix(y_pred_factored, y_test_factored)
  
  return(list(model = model, history = history, confusion_matrix = cm))
}

# ==============================================================================
#   PIPELINE COMPLETO
# ==============================================================================

if (!requireNamespace("hyperSpec", quietly = TRUE)) install.packages("hyperSpec")
library(hyperSpec)
library(caret)

# ------------------------------------------------------------------------------
# 1. CARREGAMENTO E FILTRAGEM ANTI-NA DAS CLASSES
# ------------------------------------------------------------------------------

# Opção A: Usando o pacote hyperSpec (Padrão)
if (!exists("chondro")) chondro <- hyperSpec::chondro
indices_validos <- which(!is.na(chondro$clusters))

df_chondro <- list(
  data = chondro$spc[indices_validos, ],                  
  wavenumbers = wl(chondro),            
  y_train = as.factor(chondro$clusters[indices_validos]) 
)

# [NOTA]: Se fores usar os teus ficheiros .npy locais, descomenta as linhas abaixo:
# df_bruto <- dataset_npy(x_file="X_2019clinical.npy", y_file="y_2019clinical.npy", w_file="wavenumbers.npy")
# indices_validos <- which(!is.na(df_bruto$y_train))
# df_chondro <- list(data=df_bruto$data[indices_validos, , drop=FALSE], wavenumbers=df_bruto$wavenumbers, y_train=as.factor(df_bruto$y_train[indices_validos]))


cat("Amostras limpas enviadas para o pipeline:", nrow(df_chondro$data), "\n\n")


# ------------------------------------------------------------------------------
# 2. FLUXO SEQUENCIAL DE PRÉ-PROCESSAMENTO 
# ------------------------------------------------------------------------------
cat("--- A EXECUTAR PRÉ-PROCESSAMENTO ---\n")

alinhamento_teste <- align_silicon_peak(
  wavenumber = df_chondro$wavenumbers, 
  intensity = df_chondro$data[1, ], 
  target_peak = 1003, 
  search_window = 15
)

df_cropped <- crop_espetro(df_chondro, min_wn = 650, max_wn = 1750)
df_despiked <- dataset_despike(df_cropped, ma = 10, threshold = 3.5)
df_baseline <- airPLS_fast_dataset(df_despiked, lambda = 1e4)

df_minmax <- raman_norm(df_baseline, method = "min-max")
df_l2     <- raman_norm(df_baseline, method = "L2")


# ------------------------------------------------------------------------------
# 3. EXTRAÇÃO E TRANSFORMAÇÕES MATEMÁTICAS ADICIONAIS
# ------------------------------------------------------------------------------
cat("\n--- A EXECUTAR TRANSFORMAÇÕES MATEMÁTICAS ---\n")

df_wavelet <- raman_transform_wavelet(df_l2, level = 1)
cat("Dimensões após Transformada Wavelet Haar:", dim(df_wavelet$data), "\n")

df_fourier <- raman_transform_fourier(df_l2)
cat("Dimensões após Transformada de Fourier:", dim(df_fourier$data), "\n")


# ------------------------------------------------------------------------------
# 4. BLOCO DE VISUALIZAÇÃO GRÁFICA 
# ------------------------------------------------------------------------------
cat("\n--- A GERAR PROJEÇÕES GRÁFICAS ---\n")

par(mfrow = c(2, 1))
plot(df_minmax$wavenumbers, colMeans(df_minmax$data), type="l", col="blue", lwd=1.5,
     main="Tua Normalização Min-Max (Média Chondro)", xlab="cm-1", ylab="Intensidade")
plot(df_l2$wavenumbers, colMeans(df_l2$data), type="l", col="darkgreen", lwd=1.5,
     main="Tua Normalização Vetorial L2 (Média Chondro)", xlab="cm-1", ylab="Intensidade")
par(mfrow = c(1, 1))

df_radar_indiv <- raman_radar_plot(df_baseline, sample_idx = 1)
raman_radar_class_averages(df_l2)
meu_pca_res <- raman_pca_plot(df_l2)


# ------------------------------------------------------------------------------
# 5. MOCK EXTRA: EXTRAÇÃO DE FEATURES DE PICOS
# ------------------------------------------------------------------------------
cat("\n--- A EXTRAIR FEATURES ESTRUTURAIS DE PICOS ---\n")
nomes_falsos_espectros <- paste0(df_l2$y_train, "$sample_", 1:nrow(df_l2$data))
picos_alvo <- c(1003, 1450, 1660) 

lista_propriedades_mock <- lapply(1:nrow(df_l2$data), function(i) {
  list(
    Peaks = picos_alvo,
    peak_heights = c(df_l2$data[i, which.min(abs(df_l2$wavenumbers - 1003))], 
                     df_l2$data[i, which.min(abs(df_l2$wavenumbers - 1450))], 
                     df_l2$data[i, which.min(abs(df_l2$wavenumbers - 1660))]),
    prominences = c(0.5, 0.8, 0.6), widths = c(10, 15, 12), width_heights = c(5, 7, 6),
    left_minima = c(990, 1430, 1640), right_minima = c(1015, 1470, 1680), AUC = c(5, 12, 8)
  )
})

df_features_extraidas <- create_features_dataset(
  d_spectra_names = nomes_falsos_espectros, 
  properties_list = lista_propriedades_mock, 
  vector = picos_alvo
)
cat("Matriz estatística de features criada com sucesso! Dimensões:", dim(df_features_extraidas), "\n\n")


# ------------------------------------------------------------------------------
# 6. CAPTURA E IMPRESSÃO DE RESULTADOS DE TODOS OS MODELOS PREDICTIVOS
# ------------------------------------------------------------------------------
cat("==================================================================\n")
cat("   TESTE DE MODELOS PREDICTIVOS (MACHINE LEARNING & DEEP LEARNING)\n")
cat("==================================================================\n")

# --- MODELO A: RANDOM FOREST  ---
cat("\n>>> Executando: Random Forest (Normalização L2) ...\n")
rf_l2_log <- capture.output({
  rf_l2_vis <- withVisible(ML_raman(df_l2, method = "rf", train_prop = 0.7))
})

cat("\n[RESULTADOS - RANDOM FOREST L2]:\n")
if (!is.null(rf_l2_vis$value)) {
  print(rf_l2_vis$value$results[, c("mtry", "Accuracy", "Kappa")])
  cat("\nMatriz de Confusão Acumulada nos Folds:\n")
  print(rf_l2_vis$value$resampledCM)
} else {
  cat(paste(tail(rf_l2_log, 25), collapse = "\n"), "\n")
}


# --- MODELO B: SVM RADIAL  ---
cat("\n>>> Executando: Support Vector Machine (Normalização L2) ...\n")
svm_l2_log <- capture.output({
  svm_l2_vis <- withVisible(ML_raman(df_l2, method = "svm", train_prop = 0.7))
})

cat("\n[RESULTADOS - SVM RADIAL L2]:\n")
if (!is.null(svm_l2_vis$value)) {
  print(svm_l2_vis$value$results[, c("sigma", "C", "Accuracy", "Kappa")])
  cat("\nMatriz de Confusão Acumulada nos Folds:\n")
  print(svm_l2_vis$value$resampledCM)
  
  importancia <- varImp(svm_l2_vis$value)
  plot(importancia, top = 20, main = "Top 20 Wavenumbers Críticos (SVM - L2)")
} else {
  cat(paste(tail(svm_l2_log, 25), collapse = "\n"), "\n")
}


# --- MODELO C: RANDOM FOREST ---
cat("\n>>> Executando: Random Forest (Normalização Min-Max) ...\n")
rf_mm_log <- capture.output({
  rf_mm_vis <- withVisible(ML_raman(df_minmax, method = "rf", train_prop = 0.7))
})

cat("\n[RESULTADOS - RANDOM FOREST MIN-MAX]:\n")
if (!is.null(rf_mm_vis$value)) {
  print(rf_mm_vis$value$results[, c("mtry", "Accuracy", "Kappa")])
} else {
  cat(paste(tail(rf_mm_log, 20), collapse = "\n"), "\n")
}


# --- MODELO D: SVM RADIAL  ---
cat("\n>>> Executando: Support Vector Machine (Normalização Min-Max) ...\n")
svm_mm_log <- capture.output({
  svm_mm_vis <- withVisible(ML_raman(df_minmax, method = "svm", train_prop = 0.7))
})

cat("\n[RESULTADOS - SVM RADIAL MIN-MAX]:\n")
if (!is.null(svm_mm_vis$value)) {
  print(svm_mm_vis$value$results[, c("sigma", "C", "Accuracy", "Kappa")])
} else {
  cat(paste(tail(svm_mm_log, 20), collapse = "\n"), "\n")
}


# --- MODELO E: DEEP LEARNING (CNN 1D INTERNA COM KERAS) ---
cat("\n>>> Executando: Rede Neuronal Convolucional 1D (CNN 1D) ...\n")

resultado_cnn <- ML_raman_cnn1d(
  dataset = df_l2, 
  train_prop = 0.7, 
  epochs = 50,       
  batch_size = 16
)

cat("\n[RESULTADOS - CNN 1D DEEP LEARNING]:\n")
cat("Accuracy final obtida nos dados de teste:", round(resultado_cnn$metrics * 100, 2), "%\n\n")

cat("Matriz de Confusão Detalhada (Conjunto de Validação Independente):\n")
print(resultado_cnn$confusion_matrix)

# Plota as curvas corrigidas de loss/accuracy 
plot(resultado_cnn$history)


# --- TESTAR VISUALIZAÇÃO ---
# Compara o df_chondro original com o df totalmente tratado
visualizar_raw_vs_pre(df_raw = df_chondro, df_processado = df_l2, metodo_norm = "L2")



# --- TESTAR CNN 2D COM RADAR PLOTS ---
# Executa a classificação baseada na geometria circular do radar plot
resultado_radar_cnn <- ML_raman_cnn2d_radar(
  dataset = df_l2, 
  img_size = 64,    # Resolução da imagem (64x64 pixéis)
  epochs = 20, 
  batch_size = 32
)

# Ver identificação e performance do modelo
print(resultado_radar_cnn$confusion_matrix)
plot(resultado_radar_cnn$history)
                                     
