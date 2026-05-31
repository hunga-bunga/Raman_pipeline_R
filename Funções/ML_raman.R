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
    
    # Validação cruzada que também respeita o agrupamento de individuos nos folds
    fitControl <- trainControl(method = "cv", number = 5,
                               index = groupKFold(dataset$patient_ids[trainIndex], k = 5))
  } else {
    warning("Aviso: 'patient_ids' não encontrado. A usar divisão por espetro (Risco de Data Leakage).")
    trainIndex <- createDataPartition(df_ml$class, p = train_prop, list = FALSE)
    
    train_data <- df_ml[trainIndex, ]
    test_data  <- df_ml[-trainIndex, ]
    
    fitControl <- trainControl(method = "cv", number = 5)
  }
  
  #CÁLCULO DOS PESOS PARA O BALANCEAMENTO
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