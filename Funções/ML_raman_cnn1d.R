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
