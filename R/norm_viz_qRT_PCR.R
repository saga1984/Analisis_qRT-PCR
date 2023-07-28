#' @title Obtiene Tablas y Box Plot de CTs de genes bajo diferentes tratamientos/condiciones
#'
#' @description Toma archivos de corridas de qRT-PCR. Filtra y transforma datos para obtener CTs y CTs sin valores extremos
#' 		devuelve tablas e imagenes de calidad de publicacion

#'
#' @param ruta Ruta de carpeta origen en winsys
#' @param archivos Vector de archivos de los diferentes Tratamientos y Condiciones  
#' @param vector_normalizadores vector de nombre de genes a evaluar (no necesiamente tienen que ser normalizadores)
#' @return Obtiene y filtra datos de corridas para analisis de qRT PCR (separados por corrida)
#'
#' @export


# primera funcion para analisis de qRT_PCR
viz_norm_qRT_PCR <- function(ruta, 
                            archivos,
                            vector_normalizadores) {
  
  # iterar sobre archivos
  for (archivo in archivos) {
 
    # crear tibble (removiendo espacios en blanco)
    datos <- read_excel(paste(ruta,"/", archivo, sep = ""), trim_ws = TRUE)
    
    # volver data frame
    datos <- as.data.frame(datos)
    
    ############################## filtrado ########################################
    
    # filtrar columnas (por nombre)
    datos <- datos[,c("Target", "Sample", "Biological Set Name", "Content", "Cq")]
    
    # eliminar NTC´s, filtrando columna Content
    datos <- datos[datos[,"Content"] != "NTC",]
    
    # eliminar Targets no deseados, filtrando columna Target
    datos <- subset(datos, Target %in% vector_normalizadores)
    
    # crear nueva columna de ID completo
    datos$Id <- paste(datos$Target, 
                      datos$Sample,
                      datos$`Biological Set Name`)
    
    # si hay NAs (no se especifico condicion), agregarla
    if (any(grepl("NA", datos$Id))) {
      
      # definir valores
      valores <- unique(datos[, "Biological Set Name"])
      
      # remover NAs
      valor <- valores[!is.na(valores)]
      
      # reemplazar NAs del Id
      datos$Id <- str_replace(datos$Id, "NA", valor)
    }
    
    # dejar solo dos columnas
    datos <- datos[,c("Id", "Cq")]
    
    # ordenar en base a id
    datos <- datos[order(datos[,"Id"]),]
    
    # volver Cq una columna de integers
    datos$Cq <- as.double(datos$Cq)
        
    # tratamiento_condicion
    tratamiento_condicion <- str_replace(
      gsub("[A-z]+/|_corrida[12]|.xlsx", "", archivo), 
      "_", 
      " ")
    
    # agregar columna de tratamiento_condicion
    datos$Trat_Cond <- rep(tratamiento_condicion, nrow(datos))
    
    # asiganar nombre unico en el entorno padre
    assign(
      paste("norm_viz", gsub(" ", "_", tratamiento_condicion), sep = "_"),
      datos, envir = parent.frame())
    
    
      ############################ terminar loop ###############################
  }
  
  ############################## graficar ######################################
  ############################### TODOS ########################################
  
  # unir data frames
  todos_df <- ls(envir = parent.frame())
  
  # filtrar por string, norm_viz data frames
  todos_norm_viz <- todos_df[grepl("norm_viz_", todos_df)]
  
  # crear lista de objetos
  lista_norm_viz <- mget(todos_norm_viz, envir = parent.frame())
  
  # volver data frame
  df_norm_viz <<- as.data.frame(do.call(rbind, lista_norm_viz))
  
  # determinar Cq maximo
  cq_maximo <- max(df_norm_viz$Cq)
  
  ########################### GRAFICAR TODOS ###################################
  
  # definir formatos
  # formatos <- c("tiff", "jpeg")
  formatos <- "jpeg"
  
  # si existe ruta
  if (dir.exists(ruta)) {
    # guardar imagen en formatos pre-establecidos
    for(i in formatos) {
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ruta, "/", "plot_", paste(vector_normalizadores, collapse = "_"), ".", i, sep = ""),
                   res = 300,
                   width = 8000,
                   height = 4000)
      # crear y guardar el heatmpat euclidean
      print(
        ggplot(data = df_norm_viz,
               aes(x = Id,
                   y = Cq,
                   fill = Trat_Cond)) +
          geom_boxplot(outlier.colour = "red", varwidth = TRUE) +
          geom_jitter(shape = 16, size = 0.2) +
          theme_minimal() +
          theme(legend.position = "left",
                axis.text.x = element_text(angle = 90)) +
          scale_fill_brewer(palette = "Paired") +
          ylim(0, cq_maximo + 5) +
          labs(title = "Normalization/Endogenos Genes on Abiotic Stress",
               subtitle = paste(vector_normalizadores, collapse = " & "),
               caption = "Image 5 ",
               x = "Treatment Condition", y = "Cycling Threshold", 
               fill = "")
      )
      dev.off()
    }
  }  
  
  ###################### remover valores extremos ##############################
  
  # volver Cq una columna de integers
  df_norm_viz$Cq <- as.double(df_norm_viz$Cq)
  
  # obtener tabla de promedios de Cq (o Ct)
  promedios <- aggregate(Cq ~ Id, 
                         data = df_norm_viz, 
                         FUN = mean)  
  
  # agregar promedios a cada tipo de fila
  df_norm_viz_modif <- merge(df_norm_viz, promedios, by = "Id", all.x = TRUE)
  
  # modificar columnas
  colnames(df_norm_viz_modif) <- c("Id", "Cq", "Trat_Cond", "Mean")
  
  # agregar columna de error absoluto
  df_norm_viz_modif$err_abs <- abs(df_norm_viz_modif$Mean - df_norm_viz_modif$Cq)
  
  # obtener outliers
  quitar <- aggregate(err_abs ~ Id, 
                      data = df_norm_viz_modif, 
                      FUN = max) 
  
  # quitar outliers
  df_norm_viz_modif_func <- anti_join(df_norm_viz_modif, 
                                 quitar, 
                                 by = "err_abs")
  
  
  # quitar outliers
  df_norm_viz_modif <<- df_norm_viz_modif_func
  
  # determinar Cq maximo
  cq_maximo_modif <- max(df_norm_viz_modif_func$Cq)
  
  ################## GRAFICAR TODOS QUITANDO VALORES EXTREMOS ##################
  
  # definir formatos
  # formatos <- c("tiff", "jpeg")
  formatos <- "jpeg"
  
  # si existe ruta
  if (dir.exists(ruta)) {
    # guardar imagen en formatos pre-establecidos
    for(i in formatos) {
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(ruta, "/", "plot_sin_extremos_", paste(vector_normalizadores, collapse = "_"), ".", i, sep = ""),
                   res = 300,
                   width = 8000,
                   height = 4000)
      # crear y guardar el heatmpat euclidean
      print(
        ggplot(data = df_norm_viz_modif_func,
               aes(x = Id,
                   y = Cq,
                   fill = Trat_Cond)) +
          geom_boxplot(outlier.colour = "red", varwidth = TRUE) +
          geom_jitter(shape = 16, size = 0.2) +
          theme_minimal() +
          theme(legend.position = "left",
                axis.text.x = element_text(angle = 90)) +
          scale_fill_brewer(palette = "Paired") +
          ylim(0, (cq_maximo_modif + 5)) +
          labs(title = "Normalization/Endogenos Genes on Abiotic Stress",
               subtitle = paste(vector_normalizadores, collapse = " & "),
               caption = "Image 5 ",
               x = "Treatment Condition", y = "Cycling Threshold", 
               fill = "")
      )
      dev.off()
    }
  }
  
  ###################### POR TRATAMIENTO Y POR GEN #############################
  
  #################### crear columna de solo tratamientos ######################
  
  # crear columna de solo tratamientos
  df_norm_viz_trat_cond_func <- df_norm_viz %>% 
    mutate(Trat_Cond = sapply(strsplit(Trat_Cond, " "), "[", 1))
  
  # volver global
  df_norm_viz_trat_cond <<- df_norm_viz_trat_cond_func
  
  
  ############## crear columna de solo tratamientos y genes ####################
  
  # crear columna de solo tratamientos
  df_norm_viz_gen_trat_func <- df_norm_viz_trat_cond_func %>% 
    mutate(Id = sapply(strsplit(Id, " "), "[", 1))
  
  # volver global 
  df_norm_viz_gen_trat <<- df_norm_viz_gen_trat_func
  
  ############## crear columna de solo tratamientos y genes_ctl ####################

  ########################### GRAFICAR  ###################################

  # definir formatos
  # formatos <- c("tiff", "jpeg")
  formatos <- "jpeg"

  # definir vector de archivos
  nombres <- list(df_norm_viz_trat_cond_func,
               df_norm_viz_gen_trat_func)

  # si existe ruta
  if (dir.exists(ruta)) {
    # guardar imagen en formatos pre-establecidos
    for(i in formatos) {
      for (nombre in 1:length(nombres)) {
        # crear y guardar los heatmaps
        match.fun(i)(paste(ruta, "/", "plot_", nombre, 
                           "_",
                           paste(vector_normalizadores, collapse = "_"), 
                           ".", 
                           i, 
                           sep = ""),
                     res = 300,
                     width = 8000,
                     height = 4000)
        # crear y guardar el heatmpat euclidean
        print(
          ggplot(data = nombres[[nombre]],
                 aes(x = Id,
                     y = Cq,
                     fill = Trat_Cond)) +
            geom_boxplot(outlier.colour = "red", varwidth = TRUE) +
            geom_jitter(shape = 16, size = 0.4) +
            theme_minimal() +
            theme(legend.position = "top",
                  axis.text.x = element_text(angle = 90)) +
            scale_fill_brewer(palette = "Dark2") +
            ylim(0, cq_maximo + 5) +
            labs(title = "Normalization/Endogenos Genes on Abiotic Stress",
                 subtitle = paste(vector_normalizadores, collapse = " & "),
                 caption = "Image 5 ",
                 x = "Treatment/Condition/Gene", y = "Cycling Threshold",
                 fill = "")
        )
        dev.off()
      }
    }
  }

  ########################### sin extremos #####################################
  #################### crear columna de solo tratamientos ######################
  
  # crear columna de solo tratamientos
                                      
  df_norm_viz_trat_cond_modif_func <- df_norm_viz_modif_func %>% 
    mutate(Trat_Cond = sapply(strsplit(Trat_Cond, " "), "[", 1))
  
  # volver global
  df_norm_viz_trat_cond_modif <<- df_norm_viz_trat_cond_modif_func
                                  
  

  
  ############# crear columna de solo tratamientos y genes #####################
  
  # crear columna de solo tratamientos
  df_norm_viz_gen_trat_modif_func <- df_norm_viz_trat_cond_modif_func %>% 
    mutate(Id = sapply(strsplit(Id, " "), "[", 1))
  
  # volver global 
  df_norm_viz_gen_trat_modif <<- df_norm_viz_gen_trat_modif_func
 
  ######################## graficar ##########################
  
  # formatos
  # formatos <- c("jpeg", "tiff")
  formatos <- "jpeg"
  
  # definir vector de archivos
  nombres2 <- list(df_norm_viz_trat_cond_modif_func,
                   df_norm_viz_gen_trat_modif_func)
              
          
  # si existe ruta
  if (dir.exists(ruta)) {
    # guardar imagen en formatos pre-establecidos
    for(i in formatos) {
      for (nombre_modif in 1:length(nombres2)) {
        # crear y guardar los heatmaps
        match.fun(i)(paste(ruta, "/", "plot_sin_extremos_", nombre_modif,
                           "_",
                           paste(vector_normalizadores, collapse = "_"),
                           ".", 
                           i, 
                           sep = ""),
                     res = 300,
                     width = 8000,
                     height = 4000)
        # crear y guardar el heatmpat euclidean
        print(
          ggplot(data = nombres2[[nombre_modif]],
                 aes(x = Id,
                     y = Cq,
                     fill = Trat_Cond)) +
            geom_boxplot(outlier.colour = "red", varwidth = TRUE) +
            geom_jitter(shape = 16, size = 0.4) +
            theme_minimal() +
            theme(legend.position = "top",
                  axis.text.x = element_text(angle = 90)) +
            scale_fill_brewer(palette = "Dark2") +
            ylim(0, cq_maximo + 5) +
            labs(title = "Normalization/Endogenos Genes on Abiotic Stress",
                 subtitle = paste(vector_normalizadores, collapse = " & "),
                 caption = "Image 5 ",
                 x = "Treatment/Condition/Gene", y = "Cycling Threshold",
                 fill = "")
        )
        dev.off()
      }
    }
  }

  ########################## guardar tablas #################################
  
  
  
  # unir data frames
  todos_df <- ls(envir = parent.frame())
  
  # filtrar
  tablas_norm_viz <- todos_df[grepl("df_norm_viz",  todos_df)]
  
  # Guardar cada data frame en un archivo CSV
  for (df_name in tablas_norm_viz) {
    df <- get(df_name, envir = parent.frame())
    write.csv(df, 
              file = paste(ruta, 
                           "/",
                           df_name, 
                           "_", 
                           paste(vector_normalizadores, collapse = "_"),
                           ".csv", 
                           sep = ""),
              row.names = FALSE)
    
    # imprimir que se hizo
    cat(paste("Se guardo el archivo:\n", df_name, "\n en: ", ruta, sep = ""))
    
  }
  
  
  ####################### remover data frames ##################################
  # seleccionar todos los data frames
  todos_datos <- ls(envir = parent.frame())
  
  # eliminar
  rm(list = todos_datos, envir = parent.frame())
}
