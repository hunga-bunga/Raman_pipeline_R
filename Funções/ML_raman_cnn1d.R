ML_raman_cnn1d <- function(dataset, train_prop = 0.7, epochs = 30, batch_size = 16) {
  
  # --- CARREGAMENTO AUTOMÁTICO DAS LIBRARIES ---
  # Verifica e carrega o caret
  if (!requireNamespace("caret", quietly = TRUE)) {
    message("O package 'caret' não foi encontrado. A instalar...")
    install.packages("caret")
  }
  library(caret)
  
  # Verifica e carrega o keras3 (ou keras tradicional)
  if (!requireNamespace("keras3", quietly = TRUE)) {
    message("O package 'keras3' não foi encontrado. A instalar...")
    install.packages("keras3")
    # Nota: O utilizador ainda precisará de ter o ambiente virtual do Python 
    # configurado, mas isto garante que a library de R fica pronta.
  }
  library(keras3)
  
  # ---------------------------------------------------------------------------
  
  message("A preparar os dados para a CNN 1D...")
  
  # 1. Preparação das variáveis
  X <- dataset$data
  y <- as.factor(dataset$y_train)
  
  # Codificar as labels para numérico (0, 1, 2...)
  num_classes <- length(levels(y))
  y_encoded <- as.numeric(y) - 1
  
  # 2. Divisão Treino / Teste
  set.seed(123)
  train_idx <- createDataPartition(y, p = train_prop, list = FALSE)
  
  X_train <- X[train_idx, ]
  y_train <- y_encoded[train_idx]
  X_test  <- X[-train_idx, ]
  y_test  <- y_encoded[-train_idx]
  
  # O Keras espera tensores 3D para convoluções 1D: [Amostras, Comprimento do Sinal, Canais]
  dim(X_train) <- c(nrow(X_train), ncol(X_train), 1)
  dim(X_test)  <- c(nrow(X_test), ncol(X_test), 1)
  
  # Conversão das labels para matrizes binárias (One-Hot Encoding)
  y_train_cat <- to_categorical(y_train, num_classes)
  y_test_cat  <- to_categorical(y_test, num_classes)
  
  # 3. Construção da Arquitetura da CNN 1D
  message("A construir o modelo de Deep Learning...")
  
  input_shape <- c(ncol(X), 1)
  
  model <- keras_model_sequential() %>%
    # Primeira Camada Convolucional
    layer_conv_1d(filters = 16, kernel_size = 5, activation = 'relu', input_shape = input_shape) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_dropout(rate = 0.25) %>%
    
    # Segunda Camada Convolucional
    layer_conv_1d(filters = 32, kernel_size = 3, activation = 'relu') %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_dropout(rate = 0.25) %>%
    
    # Aplanar o tensor para entrar nas camadas densas
    layer_flatten() %>%
    layer_dense(units = 64, activation = 'relu') %>%
    layer_dropout(rate = 0.5) %>%
    
    # Camada de Saída (Softmax)
    layer_dense(units = num_classes, activation = 'softmax')
  
  # Compilação do Modelo
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_adam(learning_rate = 0.001),
    metrics = c('accuracy')
  )
  
  # 4. Treino do Modelo
  message("A treinar a CNN 1D (isto pode demorar um pouco)...")
  history <- model %>% fit(
    X_train, y_train_cat,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = 0.1, 
    verbose = 1
  )
  
  # 5. Avaliação no Conjunto de Teste
  message("A avaliar o modelo nos dados de teste...")
  eval_res <- model %>% evaluate(X_test, y_test_cat, verbose = 0)
  
  # Predições para a Matriz de Confusão
  predictions_prob <- model %>% predict(X_test)
  predictions_class <- apply(predictions_prob, 1, which.max) - 1
  
  # Converter de volta para os fatores originais
  predicted_factors <- factor(levels(y)[predictions_class + 1], levels = levels(y))
  actual_factors    <- factor(levels(y)[y_test + 1], levels = levels(y))
  
  cm <- confusionMatrix(predicted_factors, actual_factors)
  
  return(list(
    model = model,
    history = history,
    test_accuracy = eval_res[["accuracy"]],
    confusion_matrix = cm
  ))
}