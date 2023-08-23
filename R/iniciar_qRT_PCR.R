#' @title Filtra y Transforma Datos de las Diferentes Corridas. Obtiene Delta Delta CT y Desviaciones Estandar (ds)
#'
#' @description Filtra y transforma datos por corrida

#'
#' @param ruta Ruta de carpeta origen en winsys
#' @param archivos Vector de archivos de las diferentes 
#' @param target Gen problema (por ejemplo cornichon CN1)
#' @param normalizador Gen normalizador (por ejemplo tubulina, Tb o factor de elongacion, Ef)
#' @param tratamiento Tratamiento (por ejemplo estres abiotico, NaCl)
#' @param tratamiento_condicion Tratamiento y condicion (por ejemplo: Tratamiento: NaCl. Condicion: 3h)
#' @return Obtiene y filtra datos de corridas para analisis de qRT PCR (separados por corrida)
#'
#' @export

# primera funcion para analisis de qRT_PCR
iniciar_qRT_PCR <- function(ruta, 
                            archivos,
                            tratamiento,
                            tratamiento_condicion,
                            target, 
                            normalizador) {
  
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
    
    # definir target y normalizador a conservar
    target_normalizador <- c(target, normalizador)
    
    # eliminar Targets no deseados, filtrando columna Target
    datos <- subset(datos, Target %in% target_normalizador)
    
    # crear nueva columna de ID completo
    datos$Id <- paste(datos$Target, 
                      datos$Sample,
                      datos$`Biological Set Name`)
    
    # remover NAs del Id
    datos$Id <- str_remove(datos$Id, "NA")
    
    # dejar solo dos columnas
    datos <- datos[,c("Id", "Cq")]
    
    # ordenar en base a id
    datos_tmp <- datos[order(datos[,"Id"]),]
    
    # eliminar filas raras (1hr y 1 HR)
    datos_tmp <- datos[!grepl("1hr|1 HR", datos[,"Id"], ignore.case = TRUE), ]
    
    # volver Cq una columna de integers
    datos_tmp$Cq <- as.double(datos_tmp$Cq)
    
    # obtener tabla de promedios de Cq (o Ct)
    promedios <- aggregate(Cq ~ Id, 
                           data = datos_tmp, 
                           FUN = mean)  
    
    # agregar promedios a cada tipo de fila
    datos_modificado <- merge(datos_tmp, promedios, by = "Id", all.x = TRUE)
    
    # modificar columnas
    colnames(datos_modificado) <- c("Id", "Cq", "Mean")
    
    # agregar columna de error absoluto
    datos_modificado$err_abs <- abs(datos_modificado$Mean - datos_modificado$Cq)
    
    # obtener outliers
    quitar <- aggregate(err_abs ~ Id, 
                        data = datos_modificado, 
                        FUN = max) 
    
    # quitar outliers
    datos <- anti_join(datos_modificado, 
                       quitar, 
                       by = "err_abs")
    
    # eliminar columnas innecesarias
    datos <- datos[,c("Id","Cq")]
    
    if(length(normalizador) > 1){
      ############## promediar para obtener super endogeno #######################
      
      ##### definir patrones a buscar con Ctl #####
      patron_normalizador <- paste(normalizador, "Ctl", sep = " ")
      
      # definir vectores logicos de presencia ausencia de patrones
      vector_logico_norm_1 <- grepl(paste(patron_normalizador, collapse = "|"), 
                                    datos$Id)
      
      # obtener subset de vector para promediar
      datos1 <- datos[vector_logico_norm_1,]
      
      # obtener promedios de pares
      pares1 <- mean(datos1$Cq[1:nrow(datos1) %% 2 == 0])
      # obtener promedios de impares
      impares1 <- mean(datos1$Cq[1:nrow(datos1) %% 2 != 0])
      
      # unir promedios de pares e impares
      promedios1 <- data.frame(Cq = rbind(pares1, impares1))
      
      # agregar columna de Id 
      promedios1$Id <- paste(paste(normalizador, collapse = " "), "Ctl", 
                             sep = " ")
      
      ####### definir patrones a buscar con tratamientos ######
      patron_normalizador <- paste(normalizador, tratamiento, sep = " ")
      
      # definir vectores logicos de presencia ausencia de patrones
      vector_logico_norm_2 <- grepl(paste(patron_normalizador, collapse = "|"), 
                                    datos$Id)
      
      # obtener subset de vector para promediar
      datos2 <- datos[vector_logico_norm_2,]
      
      # obtener promedios de pares
      pares2 <- mean(datos2$Cq[1:nrow(datos2) %% 2 == 0])
      # obtener promedios de impares
      impares2 <- mean(datos2$Cq[1:nrow(datos2) %% 2 != 0])
      
      # unir promedios de pares e impares
      promedios2 <- data.frame(Cq = rbind(pares2, impares2))
      
      # agregar columna de Id 
      promedios2$Id <- paste(paste(normalizador, collapse = " "), tratamiento, 
                             sep = " ")
      
      # union final
      super_endogeno <- rbind(promedios1, promedios2)
      
      # asignar nombre de filas de super endogeno
      # row.names(super_endogeno) <- super_endogeno$Id
      
      ################ filas de datos que contienen al target #################
      datos_target <- datos[grepl(target, datos$Id),]
      
      ################ unir final super endpgeno y target #####################
      datos <- rbind(datos_target, super_endogeno)
      
      # remover super endogeno
      rm(super_endogeno)
      
      # filtrar elementos a remover
      objects_remove <- ls()[grepl("pares", ls())]
      
      # remover elementos numerics
      rm(list = objects_remove)
    }
    
    ################################# calculos ###################################
    
    ############################# desviacion estandar ##############################
    
    # obtener tabla de promedios para graficar
    desvestas <- aggregate(Cq ~ Id, 
                           data = datos, 
                           FUN = sd)
    
    # nombrer columnas de desvestas
    colnames(desvestas) <- c("Id", "sd")
    
    # remover espacios en blanco # ya no es necesario
    # desvestas$Id <- trimws(desvestas$Id)
    
    ###################### sd promedio target/normalizador ######################
    
    assign(
      str_replace_all(paste(desvestas$Id[grepl("Ctl", desvestas$Id)], collapse = "/"), 
                      " ", "_"),
      sqrt((desvestas$sd[grepl("Ctl", desvestas$Id)][1])^2 + 
             (desvestas$sd[grepl("Ctl", desvestas$Id)][2])^2)
    )
    
    assign(
      str_replace_all(paste(desvestas$Id[grepl(paste(target, tratamiento, sep = " "), desvestas$Id)], 
                            desvestas$Id[grepl(paste(paste(normalizador, collapse = " "), tratamiento, sep = " "), desvestas$Id)], sep = "/"), 
                      " ", "_"),
      sqrt((desvestas$sd[grepl(paste(target, tratamiento, sep = " "), desvestas$Id)])^2 + 
             (desvestas$sd[grepl(paste(paste(normalizador, collapse = " "), tratamiento, sep = " "), desvestas$Id)])^2)
    )
    
    # guardar todos los elementos del entorno que sean numericos
    numeric_objects <- Filter(is.numeric, mget(ls()))
    
    # unirlos como data frame
    sd_conjunta <- as.data.frame(do.call(rbind, numeric_objects))
    
    # agregar columna de ids
    sd_conjunta$Id <- rownames(sd_conjunta)
    
    # nombrar colummnas
    colnames(sd_conjunta) <- c("sd", "id")
    
    # completar nombre a sd_conjunta para agregar nombre del archivo
    assign(paste("sd_conjunta", str_remove_all(archivo, ".xlsx"), sep = "_"),
           sd_conjunta, envir = globalenv())
    
    # completar nombre a datos para agregar nombre del archivo
    assign(paste("datos", str_remove_all(archivo, ".xlsx"), sep = "_"),
           datos, envir = globalenv())
    
    ######################## obtener deltadeltaCT  ###############################
    
    # crear nuevo data frame
    DDCT <- data.frame(matrix(nrow = 2, ncol = 4))
    
    # asignar nombre de columnas
    colnames(DDCT) <- c("DCT_NT", "DCT_T", "DDCT", "exp2DDCT")
    
    ##### filtrar datos para delta_NT #####
    
    # definir patrones a buscar
    patron_target <- paste(target, "Ctl", sep = " ")
    patron_normalizador <- paste(paste(normalizador, collapse = " "), 
                                 "Ctl", 
                                 sep = " ")
    
    # definir vectores logicos de presencia ausencia de patrones
    vector_logico_targ <- grepl(patron_target, datos$Id)
    vector_logico_norm <- grepl(patron_normalizador, datos$Id)
    
    # llenar datos de columna DCT_NT
    DDCT[,"DCT_NT"] <- datos$Cq[vector_logico_targ] - datos$Cq[vector_logico_norm]
    
    ##### filtrar datos para delta_T #####
    
    # definir patrones a buscar
    patron_target <- paste(target, tratamiento, sep = " ")
    patron_normalizador <- paste(paste(normalizador, collapse = " "), 
                                 tratamiento, 
                                 sep = " ")
    
    # definir vectores logicos de presencia ausencia de patrones
    vector_logico_targ <- grepl(patron_target, datos$Id)
    vector_logico_norm <- grepl(patron_normalizador, datos$Id)
    
    # llenar datos de columna DCT_T
    DDCT[, "DCT_T"] <- datos$Cq[vector_logico_targ] - datos$Cq[vector_logico_norm]
    
    # llenar datos de columna DCT_T
    DDCT[, "DDCT"] <- DDCT$DCT_T - DDCT$DCT_NT
    
    # llenar datos de columna exp2DDCT
    DDCT[, "exp2DDCT"] <- 2^(-DDCT$DDCT)
    
    # definir nombre filas
    nombreFilas <- paste(str_remove(archivo, ".xlsx"), 
                         paste(target, normalizador, sep = "_"),
                         sep = "_")
    
    # asiganr nombre filas
    DDCT$ID <- nombreFilas
    
    # nombre largo
    assign(paste("DDCT",
                 str_remove(archivo, ".xlsx"), 
                 paste(target, paste(normalizador, collapse = "_"), sep = "_"),
                 sep = "_"), 
           DDCT, envir = globalenv())
    
    ######################### crear directorio #################################
    
    # asignar subcarpeta para guardar posteriores resultados
    if(length(normalizador) > 1){
      directorio_final <- paste(
        ruta, "/", "SuperEndogeno_", tratamiento_condicion, "_", target, "_", paste(normalizador, collapse = "_"), sep = "")
    } else {
      directorio_final <- paste(
        ruta, "/", tratamiento_condicion, "_", target, "_", normalizador, sep = "")
    }
    
    # crear directorio
    dir.create(directorio_final)
    
  }
}