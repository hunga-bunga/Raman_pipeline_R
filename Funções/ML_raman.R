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