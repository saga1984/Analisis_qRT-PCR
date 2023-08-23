#' @title Hace Analisis Estadistico de ANOVA y TUKEY. Hace Grafica (de barras) de Tratamiento/Condicion vs Control/Normalizador con Desviacion Estandar Conjunta de Calidad de Publicacion
#'
#' @description Trabaja sobre objetos de R previamente creados.
#' Hace analisis estadistico de ANOVA, guarda resultados en winsys
#' Hace analisis estadistico de TUKEY, guarda plot en winsys
#' Hace GRAFICO DE BARRAS de tratamiento vs control, con desviacion estandar conjunta, guarda resultados en winsys
#' 
#' @param ruta Ruta de carpeta origen en winsys
#' @param tratamiento_condicion Tratamiento y condicion (por ejemplo: Tratamiento: NaCl. Condicion: 3h)
#' @param target Gen problema (por ejemplo cornichon CN1)
#' @param normalizador Gen normalizador (por ejemplo tubulina, Tb o factor de elongacion, Ef)
#' @param tratamiento Tratamiento (por ejemplo estres abiotico, NaCl)
#' @return Graficos y resumen de resultados de anova 
#' 
#' @export



# tercera funcion para analisis de qRT PCR
anova_tukey_graficas_qRT_PCR <- function(ruta,
                                         tratamiento_condicion,
                                         target,
                                         normalizador,
                                         tratamiento){

  ################################# ANOVA ########################################
  
  # asignar subcarpeta para guardar resultados
  if(length(normalizador) > 1){
    directorio_final <- paste(
      ruta, "/", "SuperEndogeno_", tratamiento_condicion, "_", target, "_", paste(normalizador, collapse = "_"), sep = "")
  } else {
    directorio_final <- paste(
      ruta, "/", tratamiento_condicion, "_", target, "_", normalizador, sep = "")
  }
  
  # obtener ANOVA para funcion
  anova_DDCT_combinados <- aov(exp2DDCT_promedio ~ ID, 
                               data = DDCT_combinados_finales)
  
  # guardar texto de resultado de analisis de anova
  summary_anova_DDCT_combinados <- capture.output(summary(anova_DDCT_combinados))
  
  # si subdirectorio existe
  if (dir.exists(directorio_final)) {
    # guardar ANOVA en winsys
    write.table(summary_anova_DDCT_combinados, 
                file = paste(directorio_final, 
                             "/",
                             "summary_anova_", 
                             target, 
                             "_",
                             paste(normalizador, collapse = "_"), 
                             ".txt", 
                             sep = ""), 
                sep = "\t",
                row.names = FALSE, 
                col.names = FALSE)
  }
  
  cat(paste("Se acaban de guardar lel resultado de anova de: ",
              target, "_",paste(normalizador, collapse = "_"), "_", tratamiento_condicion,  
              " en:\n", directorio_final, "\n\n",sep = ""))
  
  ################################# TUKEY ########################################
  
  # hacer prueba de tukey
  tukey_DDCT_combinados <- TukeyHSD(x = anova_DDCT_combinados, conf.level = 0.95)

  ymin <- min(anova_DDCT_combinados$effects)
  ymax <- max(anova_DDCT_combinados$coefficients)
    
  # si existe subdirectorio
  if (dir.exists(directorio_final)) {
    
    # definir formatos
    # formatos <- c("tiff", "jpeg")
    formatos <- c("jpeg")
    
    # guardar imagen en formatos pre-establecidos
    for(i in formatos) {
      
      # crear y guardar los heatmaps
      match.fun(i)(paste(directorio_final, 
                         "/",
                         "Tukey_", 
                         target, 
                         "_",
                         paste(normalizador, collapse = "_"), 
                         ".", i, sep = ""),
                   res = 600,
                   width = 10000,
                   height = 14000)
      # crear y guardar el heatmpat euclidean
      print(
        
        plot(x = tukey_DDCT_combinados,
             cex.axis=1.1,
             las = 3,
             col = "blue",
             sub = paste("tratamiento vs control", target, paste(normalizador, collapse = "_"), sep = " "),
             xlim = c((ymin - 1), (ymax + 1))
        )
      )
      #linea roja en x = 0
      abline(v = 0, col="red")
      
      dev.off()
      
    }
  }
  
  cat(paste("Se acaban de guardar la grafica de Tukey de: ",
              target, "_", paste(normalizador, collapse = "_"), "_", tratamiento_condicion,
              " en:\n", directorio_final, "\n\n",sep = ""))
  
  ##############################################################################
  ################################# graficar ###################################
  ##############################################################################
  
  # obtener tabla de promedios para graficar
  promedio_grafica <- aggregate(exp2DDCT_promedio ~ ID, 
                                data = DDCT_combinados_finales, 
                                FUN = mean)
  
  
  ################ crear columna con strings para agrupar ########################
  
  # crear columna con strings para agrupar
  grupo <- data.frame(matrix(nrow = 2, ncol = 1))
  
  # indices de tratamiento. control
  TARGET <- which(grepl(target, promedio_grafica$ID))
  CNTRL <- which(grepl("Ctl", promedio_grafica$ID))
  
  # definir indices para agrupar
  grupo[TARGET,] <- "Treatment"
  grupo[CNTRL,] <- "Control"
  
  # agregar nombre columna para agrupar
  colnames(grupo) <- "grupo"
  
  # unir data frames 
  promedio_grafica <- cbind(promedio_grafica, grupo)
  
  ###################### agregar columna de  errores #############################
  
  # crear columna con strings para agregar error
  error <- data.frame(matrix(nrow = 2, ncol = 1))
  
  # indices de tratamiento. control
  TARGET <- which(grepl(paste(target, tratamiento), sd_combinados_promedio$id))
  CNTRL <- which(grepl(paste(target, "Ctl"), sd_combinados_promedio$id))
  
  # definir indices para agrupar
  error[TARGET,] <- sd_combinados_promedio$sd[TARGET]
  error[CNTRL,] <- sd_combinados_promedio$sd[CNTRL]
  
  # agregar nombre columna para agrupar
  colnames(error) <- "error"
  
  # unir data frames 
  promedio_grafica <- cbind(promedio_grafica, error)
  
  # crear tabla final para graficar con nombre largo
  assign(paste("tabla_grafica_", 
               tratamiento_condicion, 
               "_", target, 
               "_", 
               paste(normalizador, collapse = "_"), 
               sep = ""),
         promedio_grafica, envir = globalenv())
  
  ################################## GRAFICAR ####################################
  
  # sumar error y exp2DDCT_promedio
  maximos <- promedio_grafica$exp2DDCT_promedio + promedio_grafica$error
  
  # encontrar el valor del eje y mas alto
  ymax <- max(maximos)
  
  # definir formatos
  # formatos <- c("tiff", "jpeg")
  formatos <- c("jpeg")
  
  # si existe el subdirectorio
  if (dir.exists(directorio_final)) {
    # guardar imagen en formatos pre-establecidos
    for(i in formatos) {
      
      # crear y guardar graficas de expresion relativa
      match.fun(i)(paste(directorio_final, "/", "plot_", target, "_", paste(normalizador, collapse = "_"), ".", i, sep = ""),
                   res = 300,
                   width = 5000,
                   height = 7000)
      # crear y guardar el heatmpat euclidean
      print(
        # graficar
        ggplot(data = promedio_grafica, 
               mapping = aes(
                 x = ID,
                 y = exp2DDCT_promedio, 
                 fill = grupo)) +
          geom_bar(stat = "identity", position = "stack") +
          # barras de error
          geom_errorbar(mapping = aes(ymin = exp2DDCT_promedio - error/2,
                                      ymax = exp2DDCT_promedio + error/2),
                        width=.2,
                        position=position_dodge(.9)) +
          theme_minimal() +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45),
                axis.text = element_text(size = 20)) +
          labs(x = "", y = "") +
          # cambiar colores
          scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
          ylim(0, (ymax + (0.1 * ymax)))
      )
      dev.off()
    }
    
    cat(paste("Se acaban de guardar la grafica de expresion relativa de: ",
              target, "_", paste(normalizador, collapse = "_"), "_", tratamiento_condicion, " en:\n",
              directorio_final, "\n\n", sep = ""))
  }
  
}

