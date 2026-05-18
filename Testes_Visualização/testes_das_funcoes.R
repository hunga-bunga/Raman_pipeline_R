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

# 2. Aplicar as tuas funções "de raiz" (AIRPLS e DRPLS)
# O parâmetro smoothing (lambda) controla a rigidez da linha de base.
# Valores entre 1e4 e 1e6 costumam funcionar muito bem.
valor_smoothing <- 1e5

res_airpls <- AIRPLS(wavenumber = x_espectro, intensity = y_espectro, smoothing = valor_smoothing)
res_drpls  <- DRPLS(wavenumber = x_espectro, intensity = y_espectro, smoothing = valor_smoothing)

# 3. Criar os objetos de dataset corrigidos (para manter a estrutura do teu script)
# Se quiseres avançar com um deles para a normalização, basta escolher qual usar
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

# 1. Aplicar os dois métodos ao teu dataset (partindo do df_baseline)
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

# Exemplo: Testar com Support Vector Machine (svm)
resultado_svm <- ML_raman(df_minmax, method = "svm")

# 1. Ver a Accuracy e Kappa
print(resultado_svm$results$Accuracy)

# 2. Ver a Matriz de Confusão (quais classes ele está a confundir)
# Extrai a matriz de confusão acumulada das 5 reamostragens (folds)
print(resultado_svm$resampledCM)

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

# Se quiseres ver o gráfico radar da Amostra 5, por exemplo:
df_radar <- raman_radar_plot(df_baseline, sample_idx = 5)
