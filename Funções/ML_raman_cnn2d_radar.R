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