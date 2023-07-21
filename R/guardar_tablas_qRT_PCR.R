#' @title Guarda en el Sistema de Archivos de Windows los Archivos Principales Generados en R Durante (y para) el Analisis de Expresion Relativa por qRT-PCR
#'
#' @description Guarda los objetos de R principales, en winsys
#' Incluyen:
#' Tablas de desviacion estandar de todas las corridas y promediadas
#' Tablas usadas para graficas
#' Tablas usadas para analisis estadistico de ANOVA y Tukey
#' Tablas de desviaciones estandar (ds) por corrida, combinadas (por tratamiento y contol) y ds promedio
#'
#' @param ruta Ruta de carpeta origen en winsys
#' @param tratamiento_condicion Tratamiento y condicion (por ejemplo: Tratamiento: NaCl. Condicion: 3h)
#' @param target Gen problema (por ejemplo cornichon CN1)
#' @param normalizador Gen normalizador (por ejemplo tubulina, Tb o factor de elongacion, Ef)
#' @return Graficos y resumen de resultados de anova 
#'
#' @export


# cuarta funcion para analaisis de qRT PCR
guardar_tablas_qRT_PCR <- function(ruta,
                                   tratamiento_condicion,
                                   target,
                                   normalizador){
  
  # asignar subcarpeta para guardar resultados
  directorio_final <- paste(
    ruta, "/", tratamiento_condicion, 
    "_", target, "_", normalizador,
    sep = "")
  
  if (dir.exists(directorio_final)) {
  
    # crear un vector de objetos existentes en el parent (calling) environment 
    df_names <- ls(envir = parent.frame())
    
    # eliminar objetos de anova_ de df_names
    df_names <- df_names[!grepl("anova_", df_names)]
    
    # Obtener los nombres de todos los data frames en el entorno de trabajo
    # df_names <- Filter(is.data.frame, get(df_names))
    
    #########################  filtrar DDCT (guardar todos) #####################
    
    clave_DDCT <- "DDCT_"
    
    # filtrar
    df_DDCT <- df_names[grepl(clave_DDCT, df_names)]
    
    # Guardar cada data frame en un archivo CSV
    for (df_name in df_DDCT) {
      df <- get(df_name, envir = parent.frame())
      write.csv(df, 
                file = paste(directorio_final, 
                             "/",
                             df_name, 
                             "_", 
                             tratamiento_condicion,
                             ".csv", 
                             sep = ""),
                row.names = FALSE)
    }
    
    ##################### filtrar sd_ (guardar solo combinados) #################
    
    clave_sd <- "sd_combinados"
    
    # filtrar
    df_sd <- df_names[grepl(clave_sd, df_names)]
    
    # Guardar cada data frame en un archivo CSV
    for (df_name in df_sd) {
      df <- get(df_name, envir = parent.frame())
      write.csv(df, 
                file = paste(directorio_final, 
                             "/",
                             df_name, 
                             "_", 
                             tratamiento_condicion,
                             ".csv", 
                             sep = ""),
                row.names = FALSE)
    }
    
    ################## filtrar datos_ (guardar solo combinados) ################
    
    clave_datos <- "datos_combinados"
    
    # filtrar
    df_datos <- df_names[grepl(clave_datos, df_names)]
    
    # Guardar cada data frame en un archivo CSV
    for (df_name in df_datos) {
      df <- get(df_name, envir = parent.frame())
      write.csv(df, 
                file = paste(directorio_final, 
                             "/",
                             df_name, 
                             "_", 
                             tratamiento_condicion,
                             ".csv", 
                             sep = ""),
                row.names = FALSE)
    }
    
    ######################  filtrar tabla para graficar ########################
    
    clave_grafica <- "tabla_grafica_"
    
    # filtrar
    df_grafica <-  df_names[grepl(clave_grafica, df_names)]
    
    # Guardar cada data frame en un archivo CSV
    for (df_name in df_grafica) {
      df <- get(df_name, envir = parent.frame())
      write.csv(df, 
                file = paste(directorio_final, 
                             "/",
                             df_name, 
                             "_", 
                             tratamiento_condicion,
                             ".csv", 
                             sep = ""),
                row.names = FALSE)
    }
    
    ############# remover archivos, menos algunos seleccionados) ###############
    
    # eliminar todas las funciones de df_names
    df_names <- df_names[!grepl("qRT_PCR", df_names)]
    
    # eliminar todos los elementos de df_names
    rm(list = df_names, envir = parent.frame())
    
    # imprimir que se hizo
    cat(paste("Se guardo el archivo:\n", df_names, "\n en: ", directorio_final, sep = ""))
    
  }

}
