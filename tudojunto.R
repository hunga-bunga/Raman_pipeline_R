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
  
  # Atualizar o dataset com os novos dados e retornar
  dataset$data <- norm_data
  return(dataset)
}

#PCA E GRAFICO
raman_pca_plot <- function(dataset) {
  message("A calcular o PCA e a gerar o gráfico de clusters...")
  
  # 1. Executar o PCA na matriz de dados
  # scale. = TRUE é importante se os dados não estiverem perfeitamente normalizados
  pca_res <- prcomp(dataset$data, center = TRUE, scale. = FALSE)
  
  # 2. Calcular a variância explicada pelas duas primeiras componentes (PC1 e PC2)
  var_explicada <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  lbl_pc1 <- paste0("PC1 (", round(var_explicada[1], 1), "%)")
  lbl_pc2 <- paste0("PC2 (", round(var_explicada[2], 1), "%)")
  
  # 3. Criar um data frame com as coordenadas dos pontos e as etiquetas originais
  df_pca <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Classe = as.factor(dataset$y_train)
  )
  
  # 4. Desenhar o gráfico de dispersão (Scatter Plot) nativo do R
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
  
  # 1. Preparação dos dados
  # Assumindo que:
  # dataset$data são os espectros (X)
  # dataset$y_train são as classes/etiquetas (y)
  
  df_ml <- as.data.frame(dataset$data)
  df_ml$class <- as.factor(dataset$y_train)
  
  # 2. Divisão Treino/Teste
  set.seed(123) # Para reprodutibilidade
  trainIndex <- createDataPartition(df_ml$class, p = train_prop, list = FALSE)
  train_data <- df_ml[trainIndex, ]
  test_data  <- df_ml[-trainIndex, ]
  
  # 3. Configuração do Treino (Cross-Validation)
  # Usamos 5-fold CV para validar o modelo internamente
  fitControl <- trainControl(method = "cv", number = 5)
  
  # 4. Escolha do Algoritmo
  message(paste("A treinar modelo:", method, "..."))
  
  modelo <- switch(method,
                   "svm" = train(class ~ ., data = train_data, method = "svmRadial", trControl = fitControl),
                   "rf"  = train(class ~ ., data = train_data, method = "rf", trControl = fitControl),
                   "gbm" = train(class ~ ., data = train_data, method = "gbm", trControl = fitControl, verbose = FALSE),
                   stop("Método inválido! Escolha 'svm', 'rf' ou 'gbm'.")
  )

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
  
  # --- CARREGAMENTO AUTOMÁTICO E SEGURO DAS LIBRARIES ---
  if (!requireNamespace("caret", quietly = TRUE)) {
    message("O package 'caret' não foi encontrado. A instalar...")
    install.packages("caret")
  }
  library(caret)
  
  if (!requireNamespace("keras3", quietly = TRUE)) {
    message("O package 'keras3' não foi encontrado. A instalar...")
    install.packages("keras3")
  }
  library(keras3)
  # ---------------------------------------------------------------------------
  
  message("A preparar os dados para a CNN 1D...")
  
  # 1. Preparação e Formatação dos Dados (Tensores)
  X <- dataset$data
  y <- as.numeric(as.factor(dataset$y_train)) - 1 # Converte classes para índices 0, 1, 2...
  num_classes <- length(unique(y))
  
  # Divisão Treino / Teste (Mantendo a proporção das classes)
  set.seed(123)
  trainIndex <- createDataPartition(as.factor(y), p = train_prop, list = FALSE)
  
  X_train <- X[trainIndex, , drop = FALSE]
  y_train <- y[trainIndex]
  X_test  <- X[-trainIndex, , drop = FALSE]
  y_test  <- y[-trainIndex]
  
  # Formato exigido pelas CNNs 1D no Keras: [amostras, wavenumbers, canais]
  X_train_reshaped <- array(X_train, dim = c(nrow(X_train), ncol(X_train), 1))
  X_test_reshaped  <- array(X_test, dim = c(nrow(X_test), ncol(X_test), 1))
  
  # Conversão das labels para codificação One-Hot
  y_train_categorical <- to_categorical(y_train, num_classes)
  y_test_categorical  <- to_categorical(y_test, num_classes)
  
  input_shape <- c(ncol(X_train), 1)
  
  # 2. CONSTRUÇÃO DA ARQUITETURA DA CNN 1D (CORRIGIDA CONTRA OVERFITTING)
  message("A construir o modelo Keras com técnicas de regularização...")
  
  model <- keras_model_sequential() %>%
    # Primeira Camada Convolucional + Regularização L2 para limitar pesos gigantes
    layer_conv_1d(filters = 16, kernel_size = 5, activation = 'relu', 
                  kernel_regularizer = regularizer_l2(0.001), input_shape = input_shape) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_dropout(rate = 0.4) %>% # Dropout de 40% impede que a rede decore ruído
    
    # Segunda Camada Convolucional
    layer_conv_1d(filters = 32, kernel_size = 3, activation = 'relu',
                  kernel_regularizer = regularizer_l2(0.001)) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_dropout(rate = 0.4) %>%
    
    # Aplanar a matriz para entrar na camada densa
    layer_flatten() %>%
    
    # Camada Densa Intermédia (Reduzida para limitar a capacidade de memorização)
    layer_dense(units = 32, activation = 'relu', kernel_regularizer = regularizer_l2(0.001)) %>%
    layer_dropout(rate = 0.5) %>% # Dropout estrito de 50% antes da decisão final
    
    # Camada de Saída Softmax (Classificação Multiclasse)
    layer_dense(units = num_classes, activation = 'softmax')
  
  # Compilação do modelo
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(learning_rate = 0.001),
    metrics = c('accuracy')
  )
  
  # 3. CONFIGURAÇÃO DO EARLY STOPPING (Travar o treino no ponto ótimo)
  early_stop <- callback_early_stopping(
    monitor = "val_loss", 
    patience = 3,                 # Espera no máximo 3 épocas sem melhorias antes de parar
    restore_best_weights = TRUE   # Restaura os pesos da melhor época anterior ao overfitting
  )
  
  # 4. TREINO DO MODELO
  message("A iniciar o treino do modelo...")
  history <- model %>% fit(
    X_train_reshaped, y_train_categorical,
    epochs = epochs,
    batch_size = batch_size,
    validation_data = list(X_test_reshaped, y_test_categorical),
    callbacks = list(early_stop), # Ativa a paragem antecipada aqui
    verbose = 1
  )
  
  # 5. AVALIAÇÃO FINAL NO CONJUNTO DE TESTE
  message("A gerar previsões e a calcular a Matriz de Confusão...")
  
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

##############################################################################
#teste da função de carregamento
df <- dataset_npy(
  x_file = "X_2019clinical.npy", 
  y_file = "y_2019clinical.npy", 
  w_file = "wavenumbers.npy"
)

#teste da função de crop

df_crop <- crop_espetro(df, min_wn = 500, max_wn = 1700)

# --- VISUALIZAÇÃO ---
par(mfrow = c(2, 1))

# Gráfico Original (usando os novos nomes sincronizados)
plot(df$wavenumbers, df$data[1,], type="l", 
     main="Original (Completo)", xlab="cm-1", ylab="Intensidade", col="gray")

# Gráfico Cortado
plot(df_crop$wavenumbers, df_crop$data[1,], type="l", 
     main="Resultado do Cropping (500 - 1700)", xlab="cm-1", ylab="Intensidade", col="blue")

par(mfrow = c(1, 1))

cat("Dimensões ANTES:", dim(df$data), "\n")
cat("Dimensões DEPOIS:", dim(df_crop$data), "\n")




df_clean <- dataset_despike(df_crop, ma = 10, threshold = 3.5)

cat("Mudanças totais após Despiking:", sum(df_crop$data != df_clean$data), "pontos.\n")

par(mfrow = c(2, 1))

# Comparação numa amostra específica (ex: amostra 1)
idx <- 1
plot(df_crop$wavenumbers, df_crop$data[idx,], type="l", col="gray80", 
     main=paste("Despiking: Amostra", idx), xlab="cm-1", ylab="Intensidade")
lines(df_clean$wavenumbers, df_clean$data[idx,], col="red", lwd=1)
legend("topright", legend=c("Antes", "Depois"), col=c("gray80", "red"), lty=1, bty="n")

# Verificação global (Média dos espectros)
plot(df_clean$wavenumbers, colMeans(df_clean$data), type="l", col="blue",
     main="Média Global (Dataset Processado)", xlab="cm-1", ylab="Intensidade Média")

par(mfrow = c(1, 1))


#Teste SGT (Visualizar efeito)

# 1. Aplicar o Savitzky-Golay
# Window maior (ex: 15) torna a suavização mais visível se houver muito ruído
df_sg <- savitzky_golay(df_clean, p.order = 2, window = 15, deriv = 0)

# 2. Extrair dados para facilitar o plot
amostra_idx <- 1
x <- df_clean$wavenumbers
y_original <- df_clean$data[amostra_idx, ]
y_suavizado <- df_sg$data[amostra_idx, ]

# 3. Configurar o gráfico
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

# --- GRÁFICO 1: VISÃO GERAL (Contraste Total) ---
plot(x, y_original, type = "l", col = "gray75", lwd = 1,
     main = paste("Denoising Raman: Amostra", amostra_idx),
     xlab = "Raman Shift (cm-1)", ylab = "Intensidade")

# Linha de cor forte para o resultado
lines(x, y_suavizado, col = "#FF3300", lwd = 2) # Vermelho vibrante

legend("topright", legend = c("Original (com ruído)", "SGT Processado"),
       col = c("gray75", "#FF3300"), lty = 1, lwd = c(1, 2), bty = "n")

# --- GRÁFICO 2: ZOOM NO RUÍDO (Cores Invertidas para detalhe) ---
# Seleciona uma janela pequena, ex: 1200 a 1400 cm-1
zoom_lim <- which(x >= 1200 & x <= 1400)

plot(x[zoom_lim], y_original[zoom_lim], type = "l", col = "deepskyblue", lwd = 1.5,
     main = "Zoom: Comparação de Alta Resolução",
     xlab = "cm-1", ylab = "Intensidade")

# Linha preta ou escura sobre o azul nota-se muito bem
lines(x[zoom_lim], y_suavizado[zoom_lim], col = "black", lwd = 2)

legend("topright", legend = c("Sinal Bruto", "Sinal Limpo"),
       col = c("deepskyblue", "black"), lty = 1, lwd = 2, bty = "n")

# Reset
par(mfrow = c(1, 1))

#BASELINE CORRETION ALS
# 1. Aplicar a Correção de Linha de Base (Baseline Correction)
# Usando o método "als" (Asymmetric Least Squares) mencionado na tua imagem
df_baseline <- baseline_correction(df_sg, method = "als")

# 2. Extrair dados para o gráfico
idx <- 1
x <- df_baseline$wavenumbers
y_antes <- df_sg$data[idx, ]      # Dados suavizados, mas com baseline
y_depois <- df_baseline$data[idx, ] # Dados finais (limpos)

# 3. Criar a visualização comparativa
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

# --- GRÁFICO 1: O PROCESSO DE REMOÇÃO ---
# Mostra o espectro original e a "linha" que o algoritmo identificou como fundo
plot(x, y_antes, type = "l", col = "gray70", 
     main = "Detecção de Baseline (Fluorescência)",
     xlab = "cm-1", ylab = "Intensidade")

# O fundo removido é a diferença entre o antes e o depois
baseline_estimada <- y_antes - y_depois
lines(x, baseline_estimada, col = "red", lwd = 2, lty = 2)

legend("topright", legend = c("Sinal com Fluorescência", "Baseline Estimada (ALS)"),
       col = c("gray70", "red"), lty = c(1, 2), bty = "n")

# --- GRÁFICO 2: RESULTADO FINAL ---
# Espectro "plano", pronto para análise de Machine Learning
plot(x, y_depois, type = "l", col = "black", lwd = 1.5,
     main = "Espectro Final (Baseline Corrigida)",
     xlab = "cm-1", ylab = "Intensidade Net")

abline(h = 0, col = "white", lty = 3) # Linha no zero para referência

legend("topright", legend = "Sinal Raman Puro",
       col = "white", lty = 1, bty = "n")

# Restaurar layout
par(mfrow = c(1, 1))


# TESTE DOS MÉTODOS DE BASELINE CORRECTION (AIRPLS E DRPLS DE RAIZ)

# 1. Selecionar uma amostra para o teste (ex: amostra 1)
idx_amostra <- 1
x_espectro  <- df_sg$wavenumbers
y_espectro  <- df_sg$data[idx_amostra, ]

# 2. Aplicar as funções (AIRPLS e DRPLS)
# O parâmetro smoothing (lambda) controla a rigidez da linha de base.
# Valores entre 1e4 e 1e6 costumam funcionar bem.
valor_smoothing <- 1e5

res_airpls <- AIRPLS(wavenumber = x_espectro, intensity = y_espectro, smoothing = valor_smoothing)
res_drpls  <- DRPLS(wavenumber = x_espectro, intensity = y_espectro, smoothing = valor_smoothing)

# 3. Criar os objetos de dataset corrigidos (para manter a estrutura do teu script)
df_airpls_corrected <- df_sg
df_airpls_corrected$data[idx_amostra, ] <- res_airpls$corrected$y

df_drpls_corrected <- df_sg
df_drpls_corrected$data[idx_amostra, ] <- res_drpls$corrected$y


# --- VISUALIZAÇÃO GRÁFICA ---
# Configura o ecrã para 3 gráficos (1 linha, 3 colunas) para comparar os métodos
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))

# --- GRÁFICO 1: Espectro Original vs Linhas de Base ---
plot(x_espectro, y_espectro, type = "l", col = "gray50", lwd = 1.5,
     main = paste("1. Linhas de Base (Amostra", idx_amostra, ")"), 
     xlab = "cm-1", ylab = "Intensidade")

# Desenha as baselines por cima do sinal original
lines(res_airpls$baseline$x, res_airpls$baseline$y, col = "red", lwd = 2)
lines(res_drpls$baseline$x, res_drpls$baseline$y, col = "blue", lwd = 2, lty = 2)

legend("topright", legend = c("Sinal Original", "airPLS Baseline", "drPLS Baseline"),
       col = c("gray50", "red", "blue"), lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.9)


# --- GRÁFICO 2: Resultado do airPLS ---
plot(res_airpls$corrected$x, res_airpls$corrected$y, type = "l", col = "red", lwd = 1.2,
     main = "2. Corrigido via airPLS", xlab = "cm-1", ylab = "Intensidade Net")
abline(h = 0, col = "gray70", lty = 3) # Linha de referência no zero


# --- GRÁFICO 3: Resultado do drPLS ---
plot(res_drpls$corrected$x, res_drpls$corrected$y, type = "l", col = "blue", lwd = 1.2,
     main = "3. Corrigido via drPLS", xlab = "cm-1", ylab = "Intensidade Net")
abline(h = 0, col = "gray70", lty = 3) # Linha de referência no zero

# Reset do layout de gráficos
par(mfrow = c(1, 1))


#Normalização

# 1. Aplicar os dois métodos ao dataset (partindo do df_baseline)
df_minmax <- raman_norm(df_baseline, method = "min-max")
df_l2     <- raman_norm(df_baseline, method = "L2")

# 2. Configurar o ecrã para 2 gráficos
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
idx <- 1  # Escolhemos a primeira amostra para o teste

# --- GRÁFICO 1: MIN-MAX ---
# Cores: Fundo cinza (original) vs Azul (Min-Max)
plot(df_baseline$wavenumbers, df_baseline$data[idx,], type="l", col="gray80",
     main=paste("Normalização Min-Max (Amostra", idx, ")"),
     xlab="cm-1", ylab="Intensidade (Escala 0-1)")
lines(df_minmax$wavenumbers, df_minmax$data[idx,], col="blue", lwd=1.5)

legend("topright", legend=c("Antes", "Min-Max"), col=c("gray80", "blue"), lty=1, bty="n")

# --- GRÁFICO 2: L2 (VETORIAL) ---
# Cores: Fundo cinza (original) vs Verde (L2)
plot(df_baseline$wavenumbers, df_baseline$data[idx,], type="l", col="gray80",
     main=paste("Normalização L2 / Vetorial (Amostra", idx, ")"),
     xlab="cm-1", ylab="Intensidade Relativa")
lines(df_l2$wavenumbers, df_l2$data[idx,], col="darkgreen", lwd=1.5)

legend("topright", legend=c("Antes", "L2 Norm"), col=c("gray80", "darkgreen"), lty=1, bty="n")

# Reset do layout
par(mfrow = c(1, 1))


#TESTE ML

# Exemplo: Testar com Random Forest (rf)
resultado_rf <- ML_raman(df_l2, method = "rf")
resultado_rf <- ML_raman(df_minmax, method = "rf")
# Exemplo: Testar com Support Vector Machine (svm)
resultado_svm <- ML_raman(df_l2, method = "svm")
resultado_svm <- ML_raman(df_minmax, method = "svm")
# 1. Ver a Accuracy e Kappa
print(resultado_rf$results$Accuracy)
print(resultado_svm$results$Accuracy)


# 2. Ver a Matriz de Confusão (quais classes ele está a confundir)
# Extrai a matriz de confusão acumulada das 5 reamostragens (folds)
print(resultado_svm$resampledCM)
print(resultado_rf$resampledCM)

# 3. Visualizar a importância das variáveis (quais wavenumbers importam mais)
# 1. Calcular a importância das variáveis
importancia <- varImp(resultado_svm)

# 2. Gerar o gráfico com os nomes corretos (V1, V2...)
plot(importancia, top = 20, 
     main = "Top 20 Variáveis (V) mais importantes para o SVM",
     xlab = "Importância Relativa")

#Teste do Radar plot

# Visualizar o gráfico radar da Amostra 1 (Padrão)
df_radar <- raman_radar_plot(df_baseline, sample_idx = 1)

#ver o gráfico radar da Amostra 5, por exemplo:
df_radar <- raman_radar_plot(df_baseline, sample_idx = 4)

# ==============================================================================
# TESTE DA CNN 1D (DEEP LEARNING) COM DADOS REAIS
# ==============================================================================

# teste da CNN 1D utilizando o dataset normalizado por L2 (df_l2)
# Nota: 15 épocas para o teste inicial não demorar uma eternidade na máquina.
# Se o modelo está a aprender bem,  aumentar mais tarde para 30 ou 40.

resultado_cnn <- ML_raman_cnn1d(
  dataset = df_l2, 
  train_prop = 0.7, 
  epochs = 15, 
  batch_size = 16
)

# 1. Imprimir a Exatidão (Accuracy) global alcançada no conjunto de teste independente
cat("\n=========================================\n")
cat("   RESULTADO DA CNN 1D (DEEP LEARNING)   \n")
cat("=========================================\n")
cat("Accuracy obtida nos dados de teste:", round(resultado_cnn$test_accuracy * 100, 2), "%\n\n")

# 2. Visualizar a Matriz de Confusão detalhada (Sensibilidade, Especificidade, etc.)
print(resultado_cnn$confusion_matrix)

# 3. Gerar o gráfico de evolução do treino
# Isto vai abrir uma janela no RStudio mostrando as curvas de Loss e Accuracy por Época.
plot(resultado_cnn$history)


#TESTE E VISUALIZAÇÂO DO PCA (CLUSTERING)
# Executamos o PCA para visualizar os clusters antes de treinar os modelos
meu_pca <- raman_pca_plot(df_l2)

# 1. Garantir que os dados passaram pela baseline e normalização
# df_l2 <- raman_norm(df_baseline, method = "L2")

# 2. Chama a função do Radar das Médias
raman_radar_class_averages(df_l2)
