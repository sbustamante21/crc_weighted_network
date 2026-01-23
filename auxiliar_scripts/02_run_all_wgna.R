library(dplyr)

# ==========================================
# CONFIGURACIÓN DE RUTAS
# ==========================================
input_dir  <- "counts_iteraciones"    # Carpeta donde guardaste los 10 archivos .csv
metadata_file <- "metadata_final.csv"  # Tu archivo de metadatos (el mismo para todos)
script_wgcna <- "wgcna_main_script.R"  # El nombre de tu script original
base_outdir <- "wgcna_results"         # Carpeta raíz para los resultados

# Crear la carpeta raíz si no existe
dir.create(base_outdir, showWarnings = FALSE)

# ==========================================
# AUTOMATIZACIÓN
# ==========================================

# 1. Listar todos los archivos de conteos (ej: iter_1_counts.csv, etc.)
# Usamos un patrón para evitar leer archivos que no sean los de interés
files <- list.files(path = input_dir, pattern = "\\.csv$", full.names = TRUE)

message("Se encontraron ", length(files), " matrices para procesar.")

for (file_path in files) {
    # 2. Generar un nombre único para la carpeta de salida basado en el archivo
    # Ejemplo: de 'counts_iteraciones/iter_1_rna.csv' extrae 'iter_1_rna'
    file_name <- basename(file_path)
    folder_name <- gsub(".csv", "", file_name)
    current_outdir <- file.path(base_outdir, folder_name)
    
    dir.create(current_outdir, showWarnings = FALSE, recursive = TRUE)
    
    message("\n>>> Procesando: ", folder_name)
    message(">>> Guardando en: ", current_outdir)
    
    # 3. Llamar al script de WGCNA usando Rscript
    # Pasamos los argumentos necesarios (asegúrate de que tu script original los reciba)
    # Si tu script usa variables fijas, podrías editarlas dinámicamente o pasarlas por commandArgs
    
    system(paste(
        "Rscript", script_wgcna, 
        "--expr", file_path, 
        "--traits", metadata_file, 
        "--outdir", current_outdir
    ))
}

message("\nProceso completo. Todas las redes han sido generadas en: ", base_outdir)