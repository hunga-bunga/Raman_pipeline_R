# Automated Spectral Data Pipeline (R & Machine Learning)

Este ecossistema bioinformático desenvolvido em R automatiza o processamento, transformação e classificação preditiva de dados espetrais complexos (como Espetroscopia Raman/Chondro), eliminando o viés humano e os estrangulamentos de tempo associados à análise laboratorial manual.

### Fluxo Geral do Pipeline

[Dados Brutos .npy] -> [Pré-Processamento (Clean & AirPLS)] -> [Transformação (Wavelets/Fourier)] -> [Modelos IA (ML Clássico & CNNs)]

### Dicionário Técnico de Funções

Carregamento de Ficheiros
* dataset_npy(data_dir, x_file, y_file, w_file): Realiza o parsing binário nativo de ficheiros Python NumPy (.npy) para o R, validando e extraindo metadados, intensidades, wavenumbers e as classes associadas.

Quimiometria e Pré-Processamento
* align_silicon_peak(wavenumber, intensity, target_peak, search_window): Corrige desvios e flutuações de calibração instrumental (eixo X) deslocando todo o espetro com base num pico padrão conhecido.
* crop_espetro(dataset, min_wn, max_wn): Recorta o espetro para isolar apenas a janela biológica de maior relevância quimiométrica (região fingerprint).
* dataset_despike(dataset, ma, threshold): Deteta e elimina artefactos físicos (raios cósmicos) recorrendo ao cálculo estatístico robusto do Modified Z-Score (baseado na Mediana e MAD) e suavização local.
* airPLS_fast_dataset(dataset, lambda, max_iter): Remove de forma ágil e automatizada os desvios óticos e eletrónicos através de mínimos quadrados parciais penalizados iterativos (airPLS), forçando a linha de base para zero.
* DRPLS(wavenumber, intensity, smoothing): Alternativa de correção da linha de base usando penalizações dinâmicas de segunda ordem, protegendo sinais biológicos autênticos muito baixos.
* raman_norm(dataset, method, target_peak, search_window): Normaliza as intensidades (eixo Y) suportando três abordagens matemáticas distintas: escala global (Min-Max), normalização vetorial Euclidiana (L2) ou proporcional a um pico de referência específico (Silicon).

Engenharia de Features e Transformações
* raman_transform_wavelet(dataset, level): Aplica a Transformada Discreta de Wavelet (DWT) com base Haar para compressão do sinal e filtragem adaptativa multifrequencial.
* raman_transform_fourier(dataset): Converte o sinal do domínio espacial para o domínio da frequência via FFT, extraindo o espetro de potência do sinal molecular.
* create_features_dataset(d_spectra_names, properties_list, vector): Extrai métricas geométricas tabulares diretamente de picos de interesse específicos, tais como altura, área sob a curva (AUC), proeminência e largura a meia altura.

Visualização Gráfica e Projeções
* visualizar_raw_vs_pre(df_raw, df_processado, metodo_norm): Desenha um painel comparativo vertical para validar a eficácia quimiométrica, sobrepondo o desvio total da população e a linha média central antes e depois do tratamento.
* raman_pca_plot(dataset): Executa a Análise de Componentes Principais (PCA) e projeta um gráfico bidimensional (PC1 vs PC2) para monitorizar visualmente a agregação de agrupamentos (clusters) biológicos espontâneos.
* raman_radar_plot(dataset, sample_idx): Converte a leitura espetral unidimensional linear para coordenadas polares circulares (2*pi radianos) numa prejução geométrica escura e brilhante.
* raman_radar_class_averages(dataset): Calcula e exibe o perfil geométrico médio de cada classe biológica lado-a-lado usando pequenos múltiplos polares para mapear biomarcadores visuais.

Inteligência Artificial e Classificação
* ML_raman(dataset, method, train_prop): Treina classificadores clássicos (Random Forest e SVM com Kernel RBF) utilizando divisão rigorosa orientada por Paciente (para evitar o vazamento de dados/Data Leakage) e matrizes de peso para mitigar dados desequilibrados.
* ML_raman_cnn1d(dataset, train_prop, epochs, batch_size): Instancia e otimiza uma Rede Neuronal Convolucional Unidimensional (CNN 1D) para detetar padrões espaciais contínuos no formato das ondas espetrais.
* ML_raman_cnn2d_radar(dataset, img_size, train_prop, epochs, batch_size): Abordagem de IA geométrica que rasteriza os gráficos polares (radar plots) em matrizes de pixéis de imagem e treina uma CNN 2D para classificar as amostras com base nas características visuais da forma circular gerada.

### Especificações do Ambiente

* Dependências: O pipeline necessita apenas do pacote specmine instalado no ambiente R.
* Robustez de Validação: Divisão de dados baseada em identificação de grupos e pesos de penalização por frequência de classes.
