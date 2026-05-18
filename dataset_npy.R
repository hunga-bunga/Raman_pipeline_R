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
