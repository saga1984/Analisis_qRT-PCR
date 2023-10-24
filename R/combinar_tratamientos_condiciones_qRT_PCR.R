#' @title Combina los Diferentes Tratamientos y Condiciones en un Analisis Global de qRT PCR. Genera sus Propias Tablas, Archivos de Resultados y Graficas de ANOVA, Expresion Relativa y Tukey
#'
#' @description Guarda los objetos de R principales cobinados por tratamiento y condicion, en winsys
#' Incluyen:
#' Tablas usadas para graficas
#' Tablas usadas para analisis estadistico de ANOVA y Tukey
#'
#' @param ruta Ruta de carpeta origen en winsys
#' @param tratamiento_condicion Tratamiento y condicion (por ejemplo: Tratamiento: NaCl. Condicion: 3h)
#' @param target Gen problema (por ejemplo cornichon CN1)
#' @param normalizador Gen normalizador (por ejemplo tubulina, Tb o factor de elongacion, Ef)
#' @return Graficos y resumen de resultados de anova 
#'
#' @export


# quinta funcion para analisis de qRT PCR
combinar_tratamientos_condiciones_qRT_PCR <- function(ruta_carpeta,
                                                      vector_subcarpetas,
                                                      tratamiento_condiciones,
                                                      target,
                                                      normalizador, 
                                                      formatos = "jpeg",
                                                      resolucion = 300){
  
  ##### crear subcarpeta para guardar resultadso de tratamientos condiciones #####
  
  # asignar subcarpeta para guardar resultados
  if(length(normalizador > 1)){
    directorio_tratamientos <- paste(
      ruta_carpeta, "/", "SuperEndogeno_", tratamiento_condiciones, "_", target, "_", paste(normalizador, collapse = "_"), sep = "")
  } else {
    directorio_tratamientos <- paste(
      ruta_carpeta, "/", tratamiento_condiciones, "_", target, "_", normalizador, sep = "")
  }
  
  # si el directorio de tratamientos existe
  if (!dir.exists(directorio_tratamientos)){
    dir.create(directorio_tratamientos)
  }
  
  ################# obtener objetos de tablas para graficas en R #################
  
  # crear objetos de R (tablas para graficar)
  for (subdirectorio in vector_subcarpetas){
    
    # obtener directorio final
    directorio <- paste(ruta_carpeta, "/", subdirectorio, sep = "")
    
    # Substring a buscar en el nombre de los archivos
    substring_buscado <- "tabla_grafica"
    
    # Obtener la lista de archivos en el directorio que coinciden con el substring
    tabla <- list.files(directorio, pattern = substring_buscado, full.names = TRUE)
    
    # crear objetos de R  
    assign(paste("tabla_grafica", subdirectorio, sep = "_"),
           read.csv(tabla),
           envir = parent.frame())
  }
  
  ####################### combinar tablas para graficar ########################
  
  # guardar objeto de nombres de data frames de desviaciones
  todos_objetos <- ls(envir = parent.frame())
  
  # remover funciones de qRT_PCR
  todos_objetos <- todos_objetos[!grepl("qRT_PCR", todos_objetos)]
  
  # filtrar por palabra clave "tabla_grafica"
  filtrados_tablas_graficas <- todos_objetos[grep("tabla_grafica_", todos_objetos)]
  
  # crear listas de objetos por clave
  tablas_grafica_combinados <- mget(filtrados_tablas_graficas, envir = parent.frame())
  
  # unirlos como data frame
  tablas_grafica_combinados_func <- as.data.frame(do.call(rbind, tablas_grafica_combinados))
  
  # guardar en otro objeto
  tablas_grafica_combinados <- tablas_grafica_combinados_func
  
  # volver global
  assign(paste("tablas_grafica_combinados", tratamiento_condiciones, sep = "_"),
         tablas_grafica_combinados,
         envir = parent.frame())
  
  ################################## GRAFICAR ####################################
  
  # sumar error y exp2DDCT_promedio
  maximos <- tablas_grafica_combinados_func$exp2DDCT_promedio + 
    tablas_grafica_combinados_func$error
  
  # sumar error y exp2DDCT_promedio
  minimos <- tablas_grafica_combinados_func$exp2DDCT_promedio - 
    tablas_grafica_combinados_func$error
  
  # encontrar el valor del eje y mas alto
  ymax <- max(maximos)
  
  # encontrar el valor del eje y mas bajo
  ymin <- min(minimos)
  
  # definir formatos
  
  # si existe el subdirectorio
  if (dir.exists(directorio_tratamientos)) {
    # guardar imagen en formatos pre-establecidos
    for(i in formatos) {
      
      # crear y guardar la grafica principal
      match.fun(i)(paste(directorio_tratamientos, "/", "plot_", target, "_",paste(normalizador, collapse = "_"), "_", tratamiento_condiciones, ".", i, sep = ""),
                   res = resolucion,
                   width = 7000,
                   height = 9000)
      # crear y guardar el heatmpat euclidean
      print(
        # graficar
        ggplot(data = tablas_grafica_combinados_func, 
               mapping = aes(
                 x = ID,
                 y = exp2DDCT_promedio, 
                 fill = Grupo)) +
          geom_bar(stat = "identity", position = "stack", color = "black") +
          # barras de error
          geom_errorbar(mapping = aes(ymin = exp2DDCT_promedio - error/2,
                                      ymax = exp2DDCT_promedio + error/2),
                        width=.2,
                        position=position_dodge(.9)) +
          theme_minimal() +
          theme(legend.position = "none",
                axis.text = element_text(size = 20)) +
          labs(x = "", y = "") +
          # cambiar colores
          scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
          ylim( ifelse((ymin - (0.1 * ymin)) > 0, 0, round((ymin - (0.1 * ymin)),digits = 0))  , 
                round((ymax + (0.1 * ymax)), digits = 0)
      )
      dev.off()
    }
    
    cat(paste("Se guardó la grafica de expresion relativa de: ",
                target, "_", paste(normalizador, collapse = "_"), "_", tratamiento_condiciones, " en:\n",
                directorio_tratamientos,
                sep = "", collapse = "\n"))
  }
 
  
  ################# obtener objetos de anovas en R #################
  
  # Substring a buscar en el nombre de los archivos
  substring_buscado <- "estadistica_df_global_"
  
  # crear objetos de R (tablas para graficar)
  for (subdirectorio in vector_subcarpetas){
    
    # obtener directorio final
    directorio <- paste(ruta_carpeta, "/", subdirectorio, sep = "")
    
    # Obtener la lista de archivos en el directorio que coinciden con el substring
    tabla_anova <- list.files(directorio, pattern = substring_buscado, full.names = TRUE)
    
    # crear objetos de R  
    assign(paste("tabla_anova", subdirectorio, sep = "_"),
           read.csv(tabla_anova),
           envir = parent.frame())
  }
  
  ############################## combinar anovas ###############################
  
  # guardar objeto de nombres de data frames de desviaciones
  todos_objetos <- ls(envir = parent.frame())
  
  # filtrar por palabra clave "sd_conjunta"
  filtrados_tablas_anova <- todos_objetos[grep("tabla_anova_", todos_objetos)]
  
  # crear listas de objetos por clave
  tablas_anova_combinados <- mget(filtrados_tablas_anova, envir = parent.frame())
  
  # unirlos como data frame
  tablas_anova_combinados_func <- as.data.frame(do.call(rbind, tablas_anova_combinados))
  
  # guardar en otro objeto
  tablas_anova_combinados <- tablas_anova_combinados_func
  
  # volver global
  assign(paste("tablas_anova_combinados", tratamiento_condiciones, sep = "_"),
         tablas_anova_combinados,
         envir = parent.frame())
  
  ############################## hacer anovas  #################################

  # obtener ANOVA para funcion
  anova_DDCT_combinados <- aov(exp2DDCT_promedio ~ ID,
                               data = tablas_anova_combinados_func)

  # guardar texto de resultado de analisis de anova
  summary_anova_DDCT_combinados <- capture.output(summary(anova_DDCT_combinados))

  # si subdirectorio existe
  if (dir.exists(directorio_tratamientos)) {
    # guardar ANOVA en winsys
    write.table(summary_anova_DDCT_combinados,
                file = paste(directorio_tratamientos,
                             "/",
                             "summary_anova_",
                             tratamiento_condiciones,
                             ".txt",
                             sep = ""),
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
  }

  cat(paste("Se acaban de guardar lel resultado de anova de: ",
            target, "_", paste(normalizador, collapse = "_"), "_", tratamiento_condiciones,
            " en:\n", directorio_tratamientos, sep = ""))

  ################################ hacer tukey ################################

  # hacer prueba de tukey
  tukey_DDCT_combinados <- TukeyHSD(x = anova_DDCT_combinados, conf.level = 0.95)

  ymin <- min(anova_DDCT_combinados$effects)
  ymax <- max(anova_DDCT_combinados$coefficients)

  # si existe subdirectorio
  if (dir.exists(directorio_tratamientos)) {

    # definir formatos

    # guardar imagen en formatos pre-establecidos
    for(i in formatos) {

      # crear y guardar los Tukey
      match.fun(i)(paste(directorio_tratamientos, "/", "Tukey_",
                         tratamiento_condiciones, ".", i, sep = ""),
                   res = resolucion,
                   width = 10000,
                   height = 22000)
      # crear y guardar el heatmpat euclidean
      print(

        plot(x = tukey_DDCT_combinados,
             cex.axis=1.1,
             las = 3,
             col = "blue",
             sub = paste("tratamiento vs control", target, paste(normalizador, collapse = " "), sep = " "),
             xlim = c((ymin - 1), (ymax + 1))
        )
      )
      #linea roja en x = 0
      abline(v = 0, col="red")

      dev.off()

    }
  }

  cat(paste("Se acaban de guardar lel resultado de tukey de: ",
            target, "_",paste(normalizador, collapse = "_"), "_", tratamiento_condiciones,
            " en:\n", directorio_tratamientos, sep = ""))

  ######################## guardar objetos combinados ##########################
  
  # guardar objeto de nombres de data frames de desviaciones
  todos_objetos <- ls(envir = parent.frame())
  
  # seleccionar objetos a guardar
  guardar_objetos <- todos_objetos[grepl("_combinados_", todos_objetos)]
  
  # Guardar cada data frame en un archivo CSV
  for (df_name in guardar_objetos) {
    # crear data frame
    df <- get(df_name, envir = parent.frame())
    # guardar data frame en winsys
    write.csv(df, 
              file = paste(directorio_tratamientos, "/", df_name, ".csv", sep = ""),
              row.names = FALSE)
    # registrar el guardado
    cat(paste("Se guardo el documento: ", df_name, " en\n",
              directorio_tratamientos, "\n\n",sep = ""))
  }
  
  ############# remover archivos, menos algunos seleccionados) ###############
  
  # eliminar todas las funciones de df_names
  tablas <- ls(envir = parent.frame())
  
  # eliminar todos los elementos de df_names
  rm(list = tablas, envir = parent.frame())
  
  # imprimir que se hizo
  cat(paste("Se guardo el archivo:\n", tablas, "\n en: ", directorio_tratamientos, sep = ""))
}
