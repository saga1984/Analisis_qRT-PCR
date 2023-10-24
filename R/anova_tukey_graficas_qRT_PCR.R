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
                                         tratamiento,
                                         formatos = "jpeg",
                                         resolucion = 300){

  # asignar subcarpeta para guardar resultados
  if(length(normalizador) > 1){
    directorio_final <- paste(
      ruta, "/", "SuperEndogeno_", tratamiento_condicion, "_", target, "_", paste(normalizador, collapse = "_"), sep = "")
  } else {
    directorio_final <- paste(
      ruta, "/", tratamiento_condicion, "_", target, "_", normalizador, sep = "")
  }
  
  ##############################################################################
  ################################# graficar ###################################
  ##############################################################################
  
  # obtener tabla de promedios para graficar
  promedio_grafica <- aggregate(exp2DDCT_promedio ~ ID, 
                                data = DDCT_combinados_finales, 
                                FUN = mean)
  
  
  ################ crear columna con strings para agrupar ########################
  
  # crear columna con strings para agrupar
  Grupo <- data.frame(matrix(nrow = 2, ncol = 1))
  
  # indices de tratamiento. control
  TARGET <- which(grepl(target, promedio_grafica$ID))
  CNTRL <- which(grepl("Ctl", promedio_grafica$ID))
  
  # definir indices para agrupar
  Grupo[TARGET,] <- "Treatment"
  Grupo[CNTRL,] <- "Control"
  
  # agregar nombre columna para agrupar
  colnames(Grupo) <- "Grupo"
  
  # unir data frames 
  promedio_grafica <- cbind(promedio_grafica, Grupo)
  
  ###################### agregar columna de  errores #############################
  
  # crear columna con strings para agregar error
  error <- data.frame(matrix(nrow = 2, ncol = 1))
  
  # indices de tratamiento. control
  TARGET <- which(grepl(paste(target, tratamiento, sep = "_"), sd_combinados_promedio$id))
  CNTRL <- which(grepl(paste(target, "Ctl", sep = "_"), sd_combinados_promedio$id))
  
  # definir indices para agrupar
  error[TARGET,] <- sd_combinados_promedio$sd[TARGET]
  error[CNTRL,] <- sd_combinados_promedio$sd[CNTRL]
  
  # agregar nombre columna para agrupar
  colnames(error) <- "error"
  
  # unir data frames 
  promedio_grafica <- cbind(promedio_grafica, error)

  # volver global
  promedio_grafica_global <<- promedio_grafica
  
  # crear tabla final para graficar con nombre largo
  assign(paste("tabla_grafica_", 
               tratamiento_condicion, 
               "_", target, 
               "_", 
               paste(normalizador, collapse = "_"), 
               sep = ""),
         promedio_grafica, envir = globalenv())
  
  ################################## GRAFICAR ####################################
  # promedio_grafica <- promedio_grafica_global
  # sumar error y exp2DDCT_promedio
  maximos <- promedio_grafica$exp2DDCT_promedio + promedio_grafica$error
  
  # sumar error y exp2DDCT_promedio
  minimos <- promedio_grafica$exp2DDCT_promedio - promedio_grafica$error
  
  # encontrar los valores del eje y mas alto y mas bajo
  ymax <- max(maximos)
  ymin <- min(minimos)
    
  # si existe el subdirectorio
  if (dir.exists(directorio_final)) {
    # guardar imagen en formatos pre-establecidos
    for(i in formatos) {
      
      # crear y guardar graficas de expresion relativa
      match.fun(i)(paste(directorio_final, "/", "plot_", target, "_", paste(normalizador, collapse = "_"), ".", i, sep = ""),
                   res = resolucion,
                   width = 5000,
                   height = 7000)
      # crear y guardar el heatmpat euclidean
      print(
        # graficar
        ggplot(data = promedio_grafica, 
               mapping = aes(
                 x = ID,
                 y = exp2DDCT_promedio, 
                 fill = Grupo)) +
          geom_bar(stat = "identity", color = "black") +
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
    
    cat(paste("Se acaban de guardar la grafica de expresion relativa de: ",
              target, "_", paste(normalizador, collapse = "_"), "_", tratamiento_condicion, " en:\n",
              directorio_final, "\n\n", sep = ""))
  }

  ##############################################################################
  ############### modificar data frame para estadisticos #######################
  ##############################################################################
  
  ##### preparar df #####
  
  # vector de nombre de columnas
  nombres_col <- c("ID", "Grupo", "exp2DDCT_promedio")
  
  # vector de numero de filas
  numero_fil <- nrow(promedio_grafica_global) * 3
  
  # crear nuevo df
  estadistica_df <- as.data.frame(matrix(ncol = length(nombres_col), 
                                         nrow = numero_fil)
  )
  
  # asignar nombres de columnas
  colnames(estadistica_df) <- nombres_col
  
  ################### crear nuevo df para estadisticos ########################
  
  # poblar primeras dos columnas del data frame estadistica
  for(i in 1:nrow(promedio_grafica_global)){
    for(j in 1:3){
      # crear indices para estadistica_df
      index <- (i - 1) * 3 + j
      
      # poblar primeras dos columnas de estadistica_df
      estadistica_df[index, "ID"] <- promedio_grafica_global$ID[i]
      estadistica_df[index, "Grupo"] <- promedio_grafica_global$Grupo[i]
      
    }
  }
  
  # crear nueva columna (exp_relatexp2DDCT_promedio) aparte del estadistica df
  # crear vector vacio
  expresiones_relativas <- numeric()
  
  # iterar sobre las columnas de df de graficas
  for (i in 1:nrow(promedio_grafica_global)) {
    # append cada nuevo valor al vector de expresiones_relativas
    expresiones_relativas <- c(expresiones_relativas,
                               promedio_grafica_global[,"exp2DDCT_promedio"][i] + promedio_grafica_global$error[i]/2,
                               promedio_grafica_global[,"exp2DDCT_promedio"][i],
                               promedio_grafica_global[,"exp2DDCT_promedio"][i] - promedio_grafica_global$error[i]/2
    )
  }
  
  # poblar tercera columna del nuevo df para estadisticos
  estadistica_df[,"exp2DDCT_promedio"] <- expresiones_relativas
  
  # volver global
  estadistica_df_global <<- estadistica_df
  
  ################################# ANOVA ########################################
  
  
  # obtener ANOVA para funcion
  anova_DDCT_combinados <- aov(exp2DDCT_promedio ~ ID, 
                               data = estadistica_df)
  
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
  
  ymin <- min(tukey_DDCT_combinados$ID)
  ymax <- max(tukey_DDCT_combinados$ID)
  
  # si existe subdirectorio
  if (dir.exists(directorio_final)) {
    
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
                   res = resolucion,
                   width = 10000,
                   height = 22000)
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
}

